//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise. All Rights Reserved.
// Author: morise [at] fc.ritsumei.ac.jp (Masanori Morise)
//
// F0 estimation based on DIO (Distributed Inline-filter Operation).
// Please see styleguide.txt to show special rules on names of variables
// and fnctions.
//-----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include "./dio.h"
#include "./matlabfunctions.h"
#include "./constant_numbers.h"

//-----------------------------------------------------------------------------
// struct for RawEventByDio()
// "negative" means "zero-crossing point going from positive to negative"
// "positive" means "zero-crossing point going from negative to positive"
//-----------------------------------------------------------------------------
typedef struct {
  double *negative_interval_locations;
  double *negative_intervals;
  int number_of_negatives;
  double *positive_interval_locations;
  double *positive_intervals;
  int number_of_positives;
  double *peak_interval_locations;
  double *peak_intervals;
  int number_of_peaks;
  double *dip_interval_locations;
  double *dip_intervals;
  int number_of_dips;
} ZeroCrossings;

namespace {
//-----------------------------------------------------------------------------
// GetDownsampledSignal() calculates the spectrum for estimation.
// This function carries out downsampling to speed up the estimation process
// and calculates the spectrum of the downsampled signal.
// This function is only used in the OrigianlDio().
//-----------------------------------------------------------------------------
void GetSpectrumForEstimation(double *x, int x_length, int fs, int y_length,
    int fft_size, int decimation_ratio, fft_complex *y_spectrum) {
  double *y = new double[fft_size];

  // Downsampling
  if (decimation_ratio != 1) {
    decimate(x, x_length, decimation_ratio, y);
  } else {
    for (int i = 0; i < x_length; ++i) y[i] = x[i];
  }

  // Removal of the DC component (y = y - mean value of y)
  double mean_y = 0.0;
  for (int i = 0; i < y_length; ++i) mean_y += y[i];
  mean_y /= y_length;
  for (int i = 0; i < y_length; ++i) y[i] -= mean_y;
  for (int i = y_length; i < fft_size; ++i) y[i] = 0.0;

  fft_plan forwardFFT = fft_plan_dft_r2c_1d(fft_size, y, y_spectrum,
      FFT_ESTIMATE);
  fft_execute(forwardFFT);

  fft_destroy_plan(forwardFFT);
  delete[] y;
}

//-----------------------------------------------------------------------------
// GetBestF0Contour() calculates the best f0 contour based on stabilities of
// all candidates. The F0 whose stability is minimum is selected.
// This function is only used in the OrigianlDio().
//-----------------------------------------------------------------------------
void GetBestF0Contour(int f0_length, double **f0_candidate_map,
    double **f0_stability_map, int number_of_bands, double *best_f0_contour) {
  double tmp;
  for (int i = 0; i < f0_length; ++i) {
    tmp = f0_stability_map[0][i];
    best_f0_contour[i] = f0_candidate_map[0][i];
    for (int j = 1; j < number_of_bands; ++j) {
      if (tmp > f0_stability_map[j][i]) {
        tmp = f0_stability_map[j][i];
        best_f0_contour[i] = f0_candidate_map[j][i];
      }
    }
  }
}

//-----------------------------------------------------------------------------
// EliminateUnnaturalChange() is the 1st step of the postprocessing.
// This function eliminates the unnatural change of f0 based on allowed_range.
// This function is only used in GetFinalF0Contour().
//-----------------------------------------------------------------------------
void EliminateUnnaturalChange(double *f0_before, int f0_length,
    int voice_range_minimum, double allowed_range, double *best_f0_contour,
    double *f0_after) {
  // Initialization
  for (int i = 0; i < voice_range_minimum; ++i) f0_before[i] = 0.0;
  for (int i = voice_range_minimum; i < f0_length - voice_range_minimum; ++i)
    f0_before[i] = best_f0_contour[i];
  for (int i = f0_length - voice_range_minimum; i < f0_length; ++i)
    f0_before[i] = 0.0;

  // Processing to prevent the jumping of f0
  for (int i = 0; i < voice_range_minimum; ++i) f0_after[i] = 0.0;
  for (int i = voice_range_minimum; i < f0_length; ++i)
    f0_after[i] = fabs((f0_before[i] - f0_before[i - 1]) /
    (world::kMySafeGuardMinimum + f0_before[i])) <
    allowed_range ? f0_before[i] : 0.0;
}

//-----------------------------------------------------------------------------
// EliminateSuspectedF0() is the 2nd step of the postprocessing.
// This function eliminates the suspected f0 in the anlaut and auslaut.
// This function is only used in GetFinalF0Contour().
//-----------------------------------------------------------------------------
void EliminateSuspectedF0(double *f0_before, int f0_length,
    int voice_range_minimum, double *f0_after) {
  for (int i = 0; i < f0_length; ++i) f0_after[i] = f0_before[i];

  for (int i = voice_range_minimum; i < f0_length; ++i) {
    for (int j = 1; j < voice_range_minimum; ++j) {
      if (f0_before[i - j] == 0 || f0_before[i + j] == 0) {
        f0_after[i] = 0.0;
        break;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// CountNumberOfVoicedSections() counts the number of voiced sections.
// This function is only used in GetFinalF0Contour().
//-----------------------------------------------------------------------------
void CountNumberOfVoicedSections(double *f0_after, int f0_length,
    int *positive_index, int *negative_index, int *positive_count,
    int *negative_count) {
  *positive_count = *negative_count = 0;
  for (int i = 1; i < f0_length; ++i) {
    if (f0_after[i] == 0 && f0_after[i - 1] != 0) {
      negative_index[(*negative_count)++] = i - 1;
    } else {
      if (f0_after[i - 1] == 0 && f0_after[i] != 0)
        positive_index[(*positive_count)++] = i;
    }
  }
}

//-----------------------------------------------------------------------------
// CorrectOneF0() corrects the f0[current_index] based on
// f0[current_index + sign].
// This function is only used in ForwardCorrection() and BackwardCorrection().
//-----------------------------------------------------------------------------
bool CorrectOneF0(double **f0_map, int number_of_candidates,
    double allowed_range, int current_index, int sign, double *f0_after) {
  double reference_value1 = f0_after[current_index] * 2 -
    f0_after[current_index + sign];
  double reference_value2 = f0_after[current_index];
  double minimum_error = std::min(fabs(reference_value1 -
    f0_map[0][current_index + sign]),
    fabs(reference_value2 - f0_map[0][current_index + sign]));
  double error_value;
  for (int i = 1; i < number_of_candidates; ++i) {
    error_value =
      std::min(fabs(reference_value1 - f0_map[i][current_index + sign]),
          fabs(reference_value2 - f0_map[i][current_index + sign]));
    if (error_value < minimum_error) {
      minimum_error = error_value;
      f0_after[current_index + sign] = f0_map[i][current_index + sign];
    }
  }
  if (std::min(minimum_error / (reference_value1 + world::kMySafeGuardMinimum),
      minimum_error / (reference_value2 + world::kMySafeGuardMinimum)) >
      allowed_range) {
    f0_after[current_index + sign] = 0.0;
    return false;
  }
  return true;
}

//-----------------------------------------------------------------------------
// ForwardCorrection() is the 4th step of the postprocessing.
// This function corrects the f0 candidates from backward to forward.
// This function is only used in GetFinalF0Contour().
//-----------------------------------------------------------------------------
void ForwardCorrection(double *f0_before, int f0_length, double **f0_map,
    int number_of_candidates, double allowed_range, int *positive_index,
    int *negative_index, int negative_count, double *f0_after) {
  for (int i = 0; i < f0_length; i++) f0_after[i] = f0_before[i];

  for (int i = 0; i < negative_count; ++i) {
    for (int j = negative_index[i]; j < f0_length - 1; ++j) {
      if (false == CorrectOneF0(f0_map, number_of_candidates, allowed_range,
        j, 1, f0_after)) break;
      if (i != negative_count && j == positive_index[i + 1] - 1) {
        negative_index[j] = j;
        break;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// BackwardCorrection() is the 5th step of the postprocessing.
// This function corrects the f0 candidates from forward to backward.
// This function is only used in GetFinalF0Contour().
//-----------------------------------------------------------------------------
void BackwardCorrection(double *f0_before, int f0_length, double **f0_map,
    int number_of_candidates, double allowed_range, int *positive_index,
    int *negative_index, int positive_count, double *f0_after) {
  for (int i = 0; i < f0_length; ++i) f0_after[i] = f0_before[i];

  for (int i = positive_count - 1; i >= 0; --i) {
    for (int j = positive_index[i] + 1; j > 1; --j) {
      if (false == CorrectOneF0(f0_map, number_of_candidates, allowed_range,
        j, -1, f0_after)) break;
      if (i != 0 && j == negative_index[i - 1] + 1) {
        positive_index[j] = j;
        break;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// EliminateInvalidVoicedSection() is the final step of the postprocessing.
// This function eliminates the voiced section whose the
// duration is under 50 msec.
// This function is only used in GetFinalF0Contour().
//-----------------------------------------------------------------------------
void EliminateInvalidVoicedSection(double *f0_before, int f0_length,
    int voice_range_minimum, double *f0_after) {
  for (int i = 0; i < f0_length; ++i) f0_after[i] = f0_before[i];

  int current_index;
  for (int i = 0; i < f0_length; ++i) {
    if (f0_before[i] == 0.0) continue;
    current_index = i;
    for (int j = current_index; j < f0_length; ++j)
      if (f0_before[j] == 0.0) {
        i = j;
        break;
      }
    if ((i - current_index) > voice_range_minimum) continue;
    for (int j = i; j >= current_index; --j) f0_after[j] = 0.0;
  }
}

//-----------------------------------------------------------------------------
// GetFinalF0Contour() calculates the optimal f0 contour based on all f0
// candidates. This is the processing after GetBestF0Contour().
// This function is only used in OriginalDio().
//-----------------------------------------------------------------------------
void GetFinalF0Contour(double frame_period, int number_of_candidates, int fs,
    double **f0_map, double *best_f0_contour, int f0_length,
    double *final_f0_contour) {
  // memo:
  // First and lat 50 msec are not used as the voiced section.
  int voice_range_minimum = static_cast<int>(0.5 + 50.0 / frame_period);
  // memo:
  // This is the tentative value.
  double allowed_range = 0.1 * frame_period / 5.0;

  double *f0_tmp1 = new double[f0_length];
  double *f0_tmp2 = new double[f0_length];

  EliminateUnnaturalChange(f0_tmp1, f0_length, voice_range_minimum,
      allowed_range, best_f0_contour, f0_tmp2);
  int *positive_index = new int[f0_length];
  int *negative_index = new int[f0_length];

  EliminateSuspectedF0(f0_tmp2, f0_length, voice_range_minimum, f0_tmp1);
  int positive_count, negative_count;
  CountNumberOfVoicedSections(f0_tmp1, f0_length, positive_index,
      negative_index, &positive_count, &negative_count);
  ForwardCorrection(f0_tmp1, f0_length, f0_map, number_of_candidates,
      allowed_range, positive_index, negative_index, negative_count, f0_tmp2);
  BackwardCorrection(f0_tmp2, f0_length, f0_map, number_of_candidates,
      allowed_range, positive_index, negative_index, positive_count, f0_tmp1);
  EliminateInvalidVoicedSection(f0_tmp1, f0_length, voice_range_minimum,
      final_f0_contour);

  delete[] f0_tmp1;
  delete[] f0_tmp2;
  delete[] positive_index;
  delete[] negative_index;
}

//-----------------------------------------------------------------------------
// NuttallWindow() calculates the coefficients of Nuttall window whose length
// is y_length.
//-----------------------------------------------------------------------------
void NuttallWindow(int y_length, double *y) {
  double tmp;
  for (int i = 0; i < y_length; ++i) {
    tmp  = (i + 1 - (y_length + 1) / 2.0) / (y_length + 1);
    y[i] = 0.355768 + 0.487396 * cos(2 * world::kPi * tmp) +
      0.144232 * cos(4.0 * world::kPi * tmp) +
      0.012604 * cos(6.0 * world::kPi * tmp);
  }
}

//-----------------------------------------------------------------------------
// GetFilteredSignal() calculates the signal that is the convolution of the
// input signal and low-pass filter.
// This function is only used in RawEventByDio()
//-----------------------------------------------------------------------------
void GetFilteredSignal(int half_average_length, int fft_size,
    fft_complex *x_spectrum, int x_length, double *filtered_signal) {
  double *low_pass_filter = new double[fft_size];
  for (int i = half_average_length * 2; i < fft_size; ++i)
    low_pass_filter[i] = 0.0;
  // Nuttall window is used as a low-pass filter.
  // Cutoff frequency depends on the window length.
  NuttallWindow(half_average_length * 4, low_pass_filter);

  fft_complex *low_pass_filter_spectrum = new fft_complex[fft_size];
  fft_plan forwardFFT = fft_plan_dft_r2c_1d(fft_size, low_pass_filter,
      low_pass_filter_spectrum, FFT_ESTIMATE);
  fft_execute(forwardFFT);

  // Convolution
  double tmp = x_spectrum[0][0] * low_pass_filter_spectrum[0][0] -
    x_spectrum[0][1] * low_pass_filter_spectrum[0][1];
  low_pass_filter_spectrum[0][1] =
    x_spectrum[0][0] * low_pass_filter_spectrum[0][1] +
    x_spectrum[0][1] * low_pass_filter_spectrum[0][0];
  low_pass_filter_spectrum[0][0] = tmp;
  for (int i = 1; i <= fft_size / 2; ++i) {
    tmp = x_spectrum[i][0] * low_pass_filter_spectrum[i][0] -
      x_spectrum[i][1] * low_pass_filter_spectrum[i][1];
    low_pass_filter_spectrum[i][1] =
      x_spectrum[i][0] * low_pass_filter_spectrum[i][1] +
      x_spectrum[i][1] * low_pass_filter_spectrum[i][0];
    low_pass_filter_spectrum[i][0] = tmp;
    low_pass_filter_spectrum[fft_size - i - 1][0] =
      low_pass_filter_spectrum[i][0];
    low_pass_filter_spectrum[fft_size - i - 1][1] =
      low_pass_filter_spectrum[i][1];
  }

  fft_plan inverseFFT = fft_plan_dft_c2r_1d(fft_size,
      low_pass_filter_spectrum, filtered_signal, FFT_ESTIMATE);
  fft_execute(inverseFFT);

  // Compensation of the delay.
  int index_bias = half_average_length * 2;
  for (int i = 0; i < x_length; ++i)
    filtered_signal[i] = filtered_signal[i + index_bias];

  fft_destroy_plan(inverseFFT);
  fft_destroy_plan(forwardFFT);
  delete[] low_pass_filter_spectrum;
  delete[] low_pass_filter;
}

//-----------------------------------------------------------------------------
// CheckEvent() returns 1, provided that the input value is over 1.
// This function is for RawEventByDio().
//-----------------------------------------------------------------------------
inline int CheckEvent(int x) {
  return x > 0 ? 1 : 0;
}

//-----------------------------------------------------------------------------
// ZeroCrossingEngine() calculates the zero crossing points from positive to
// negative. Thanks to Custom.Maid http://custom-made.seesaa.net/ (2012/8/19)
// This function is only used in RawEventByDio().
//-----------------------------------------------------------------------------
int ZeroCrossingEngine(double *x, int x_length, double fs,
    double *interval_locations, double *intervals) {
  int *negative_going_points = new int[x_length];

  for (int i = 0; i < x_length - 1; ++i)
    negative_going_points[i] = 0.0 < x[i] && x[i+1] <= 0.0 ? i + 1 : 0;
  negative_going_points[x_length - 1] = 0;

  int *edges = new int[x_length];
  int count = 0;
  for (int i = 0; i < x_length; ++i)
    if (negative_going_points[i] > 0)
      edges[count++] = negative_going_points[i];

  if (count < 2) {
    delete[] edges;
    delete[] negative_going_points;
    return 0;
  }

  double *fine_edges = new double[count];
  for (int i = 0; i < count; ++i)
    fine_edges[i] =
    edges[i] - x[edges[i] - 1] / (x[edges[i]] - x[edges[i] - 1]);

  for (int i = 0; i < count - 1; ++i) {
    intervals[i] = fs / (fine_edges[i + 1] - fine_edges[i]);
    interval_locations[i] = (fine_edges[i] + fine_edges[i + 1]) / 2.0 / fs;
  }

  delete[] fine_edges;
  delete[] edges;
  delete[] negative_going_points;
  return count;
}

//-----------------------------------------------------------------------------
// GetFourZeroCrossingIntervals() calculates four zero-crossing intervals.
// (1) Zero-crossing going from negative to positive.
// (2) Zero-crossing going from positive to negative.
// (3) Peak, and (4) dip. (3) and (4) are calculated from the zero-crossings of
// the differential of waveform.
//-----------------------------------------------------------------------------
void GetFourZeroCrossingIntervals(double *filtered_signal, int x_length,
    double fs, ZeroCrossings *zero_crossings) {
  const int kMiximumNumber = x_length / 4;
  zero_crossings->negative_interval_locations = new double[kMiximumNumber];
  zero_crossings->positive_interval_locations = new double[kMiximumNumber];
  zero_crossings->peak_interval_locations = new double[kMiximumNumber];
  zero_crossings->dip_interval_locations = new double[kMiximumNumber];
  zero_crossings->negative_intervals = new double[kMiximumNumber];
  zero_crossings->positive_intervals = new double[kMiximumNumber];
  zero_crossings->peak_intervals = new double[kMiximumNumber];
  zero_crossings->dip_intervals = new double[kMiximumNumber];

  zero_crossings->number_of_negatives = ZeroCrossingEngine(filtered_signal,
      x_length, fs, zero_crossings->negative_interval_locations,
      zero_crossings->negative_intervals);

  for (int i = 0; i < x_length; ++i) filtered_signal[i] = -filtered_signal[i];
  zero_crossings->number_of_positives = ZeroCrossingEngine(filtered_signal,
      x_length, fs, zero_crossings->positive_interval_locations,
      zero_crossings->positive_intervals);

  for (int i = 0; i < x_length - 1; ++i) filtered_signal[i] =
    filtered_signal[i] - filtered_signal[i + 1];
  zero_crossings->number_of_peaks = ZeroCrossingEngine(filtered_signal,
      x_length - 1, fs, zero_crossings->peak_interval_locations,
      zero_crossings->peak_intervals);

  for (int i = 0; i < x_length - 1; ++i)
    filtered_signal[i] = -filtered_signal[i];
  zero_crossings->number_of_dips = ZeroCrossingEngine(filtered_signal,
      x_length - 1, fs, zero_crossings->dip_interval_locations,
      zero_crossings->dip_intervals);
}

//-----------------------------------------------------------------------------
// GetF0CandidatesSub() calculates the f0 candidates and deviations.
// This is the sub-function of GetF0Candidates() and assumes the calculation.
//-----------------------------------------------------------------------------
void GetF0CandidatesSub(double **interpolated_f0_set, int time_axis_length,
    double f0_floor, double f0_ceil, double boundary_f0,
    double *f0_candidates, double *f0_deviations) {
  for (int i = 0; i < time_axis_length; ++i) {
    f0_candidates[i] = (interpolated_f0_set[0][i] +
      interpolated_f0_set[1][i] + interpolated_f0_set[2][i] +
      interpolated_f0_set[3][i]) / 4.0;

    f0_deviations[i]   = sqrt((
      (interpolated_f0_set[0][i] - f0_candidates[i]) *
      (interpolated_f0_set[0][i] - f0_candidates[i]) +
      (interpolated_f0_set[1][i] - f0_candidates[i]) *
      (interpolated_f0_set[1][i] - f0_candidates[i]) +
      (interpolated_f0_set[2][i] - f0_candidates[i]) *
      (interpolated_f0_set[2][i] - f0_candidates[i]) +
      (interpolated_f0_set[3][i] - f0_candidates[i]) *
      (interpolated_f0_set[3][i] - f0_candidates[i])) / 3.0);

    if (f0_candidates[i] > boundary_f0 ||
        f0_candidates[i] < boundary_f0 / 2.0 ||
        f0_candidates[i] > f0_ceil || f0_candidates[i] < f0_floor) {
      f0_candidates[i] = 0.0;
      f0_deviations[i] = world::kMaximumValue;
    }
  }
}

//-----------------------------------------------------------------------------
// GetF0Candidates() calculates the F0 candidates based on the zero-crossings.
// Calculation of F0 candidates is carried out in GetF0CandidatesSub().
//-----------------------------------------------------------------------------
void GetF0Candidates(const ZeroCrossings *zero_crossings, double boundary_f0,
  double f0_floor, double f0_ceil, double *time_axis, int time_axis_length,
  double *f0_candidates, double *f0_deviations) {
  if (0 == CheckEvent(zero_crossings->number_of_negatives - 2) *
      CheckEvent(zero_crossings->number_of_positives - 2) *
      CheckEvent(zero_crossings->number_of_peaks - 2) *
      CheckEvent(zero_crossings->number_of_dips - 2)) {
    for (int i = 0; i < time_axis_length; ++i) {
      f0_deviations[i] = world::kMaximumValue;
      f0_candidates[i] = 0.0;
    }
    return;
  }

  double *interpolated_f0_set[4];
  for (int i = 0; i < 4; ++i)
    interpolated_f0_set[i] = new double[time_axis_length];

  interp1(zero_crossings->negative_interval_locations,
      zero_crossings->negative_intervals,
      zero_crossings->number_of_negatives,
      time_axis, time_axis_length, interpolated_f0_set[0]);
  interp1(zero_crossings->positive_interval_locations,
      zero_crossings->positive_intervals,
      zero_crossings->number_of_positives,
      time_axis, time_axis_length, interpolated_f0_set[1]);
  interp1(zero_crossings->peak_interval_locations,
      zero_crossings->peak_intervals, zero_crossings->number_of_peaks,
      time_axis, time_axis_length, interpolated_f0_set[2]);
  interp1(zero_crossings->dip_interval_locations,
      zero_crossings->dip_intervals, zero_crossings->number_of_dips,
      time_axis, time_axis_length, interpolated_f0_set[3]);

  GetF0CandidatesSub(interpolated_f0_set, time_axis_length, f0_floor,
      f0_ceil, boundary_f0, f0_candidates, f0_deviations);
  for (int i = 0; i < 4; ++i) delete[] interpolated_f0_set[i];
}

//-----------------------------------------------------------------------------
// DestroyZeroCrossings() frees the memory of array in the struct
//-----------------------------------------------------------------------------
void DestroyZeroCrossings(ZeroCrossings *zero_crossings) {
  delete[] zero_crossings->negative_interval_locations;
  delete[] zero_crossings->positive_interval_locations;
  delete[] zero_crossings->peak_interval_locations;
  delete[] zero_crossings->dip_interval_locations;
  delete[] zero_crossings->negative_intervals;
  delete[] zero_crossings->positive_intervals;
  delete[] zero_crossings->peak_intervals;
  delete[] zero_crossings->dip_intervals;
}

//-----------------------------------------------------------------------------
// RawEventByDio() calculates the zero-crossings.
// This function is only used in OriginalDio().
//-----------------------------------------------------------------------------
void RawEventByDio(double boundary_f0, double fs, fft_complex *x_spectrum,
    int x_length, int fft_size, double f0_floor, double f0_ceil,
    double *time_axis, int time_axis_length, double *f0_deviations,
    double *f0_candidates) {
  double *filtered_signal = new double[fft_size];
  GetFilteredSignal(matlab_round(fs / boundary_f0 / 2.0), fft_size, x_spectrum,
    x_length, filtered_signal);

  ZeroCrossings zero_crossings = {0};
  GetFourZeroCrossingIntervals(filtered_signal, x_length, fs,
      &zero_crossings);

  GetF0Candidates(&zero_crossings, boundary_f0, f0_floor, f0_ceil,
      time_axis, time_axis_length, f0_candidates, f0_deviations);

  DestroyZeroCrossings(&zero_crossings);
  delete[] filtered_signal;
}

//-----------------------------------------------------------------------------
// GetF0CandidateAndStabilityMap() calculates all f0 candidates and
// their stabilities.
// This function is only used in the OrigianlDio().
//-----------------------------------------------------------------------------
void GetF0CandidateAndStabilityMap(double *boundary_f0_list,
    int number_of_bands, double fs_after_downsampling, int y_length,
    double *time_axis, int f0_length, fft_complex *y_spectrum,
    int fft_size, double f0_floor, double f0_ceil,
    double **f0_candidate_map, double **f0_stability_map) {
  double * f0_candidates = new double[f0_length];
  double * f0_deviations = new double[f0_length];

  // Calculation of the acoustics events (zero-crossing)
  for (int i = 0; i < number_of_bands; ++i) {
    RawEventByDio(boundary_f0_list[i], fs_after_downsampling, y_spectrum,
        y_length, fft_size, f0_floor, f0_ceil, time_axis, f0_length,
        f0_deviations, f0_candidates);
    for (int j = 0; j < f0_length; ++j) {
      // A way to avoid zero division
      f0_stability_map[i][j] = f0_deviations[j] /
        (f0_candidates[j] + world::kMySafeGuardMinimum);
      f0_candidate_map[i][j] = f0_candidates[j];
    }
  }
  delete[] f0_candidates;
  delete[] f0_deviations;
}

//-----------------------------------------------------------------------------
// OriginalDio() estimates the F0 based on Distributed Inline-filter Operation.
//-----------------------------------------------------------------------------
void OriginalDio(double *x, int x_length, int fs, double frame_period,
    double f0_floor, double f0_ceil, double channels_in_octave, int speed,
    double *time_axis, double *f0) {
  // Calculation of fundamental parameters
// Debug 2012/09/09 0.1.2
//  int number_of_bands = 1 + static_cast<int>(log(f0_ceil / f0_floor)
//    / world::kLog2 * channels_in_octave);
  int number_of_bands = 2 + static_cast<int>(log(f0_ceil / f0_floor) /
    world::kLog2 * channels_in_octave);
  double * boundary_f0_list = new double[number_of_bands];
  for (int i = 0; i < number_of_bands; ++i)
    boundary_f0_list[i] = f0_floor * pow(2.0, i / channels_in_octave);

  // normalization
  int decimation_ratio = std::max(std::min(speed, 12), 1);
  int y_length = (1 + static_cast<int>(x_length / decimation_ratio));
  int fft_size = GetSuitableFFTSize(y_length +
    (4 * static_cast<int>(1.0 + fs / boundary_f0_list[0] / 2.0)));

  // Calculation of the spectrum used for the f0 estimation
  fft_complex *y_spectrum = new fft_complex[fft_size];
  GetSpectrumForEstimation(x, x_length, fs, y_length, fft_size,
      decimation_ratio, y_spectrum);

  // f0map represents all F0 candidates. We can modify them.
  double **f0_candidate_map = new double *[number_of_bands];
  double **f0_stability_map = new double *[number_of_bands];
  int f0_length = GetSamplesForDIO(fs, x_length, frame_period);
  for (int i = 0; i < number_of_bands; ++i) {
    f0_candidate_map[i] = new double[f0_length];
    f0_stability_map[i] = new double[f0_length];
  }

  for (int i = 0; i < f0_length; ++i)
    time_axis[i] = i * frame_period / 1000.0;

  double fs_after_downsampling = static_cast<double>(fs) / decimation_ratio;
  GetF0CandidateAndStabilityMap(boundary_f0_list, number_of_bands,
      fs_after_downsampling, y_length, time_axis, f0_length, y_spectrum,
      fft_size, f0_floor, f0_ceil, f0_candidate_map, f0_stability_map);

  // Selection of the best value based on fundamental-ness.
  double *best_f0_contour = new double[f0_length];
  GetBestF0Contour(f0_length, f0_candidate_map, f0_stability_map,
      number_of_bands, best_f0_contour);

  // Postprocessing to find the best f0-contour.
  GetFinalF0Contour(frame_period, number_of_bands, fs, f0_candidate_map,
      best_f0_contour, f0_length, f0);

  delete[] best_f0_contour;
  delete[] y_spectrum;
  for (int i = 0; i < number_of_bands; ++i) {
    delete[] f0_stability_map[i];
    delete[] f0_candidate_map[i];
  }
  delete[] f0_stability_map;
  delete[] f0_candidate_map;
  delete[] boundary_f0_list;
}

}  // namespace

int GetSamplesForDIO(int fs, int x_length, double frame_period) {
  return static_cast<int>(x_length / static_cast<double>(fs) /
    (frame_period / 1000.0)) + 1;
}

void Dio2(double *x, int x_length, int fs, const DioOption option,
    double *time_axis, double *f0) {
  OriginalDio(x, x_length, fs, option.frame_period, option.f0_floor,
      option.f0_ceil, option.channels_in_octave, option.speed, time_axis, f0);
}

void Dio(double *x, int x_length, int fs, double frame_period,
    double *time_axis, double *f0) {
  const double kTargetFs = 4000.0;
  const double kF0Floor = 80.0;
  const double kF0Ceil = 640;
  const double kChannelsInOctave = 2.0;
  const int kDecimationRatio = static_cast<int>(fs / kTargetFs);

  OriginalDio(x, x_length, fs, frame_period, kF0Floor, kF0Ceil,
      kChannelsInOctave, kDecimationRatio, time_axis, f0);
}

void InitializeDioOption(DioOption *option) {
  // You can change default parameters.
  option->channels_in_octave = 2.0;
  option->f0_ceil = 640.0;
  option->f0_floor = 80.0;
  option->frame_period = 5;
  // You can use from 1 to 12.
  // Default value for 44.1 kHz of fs.
  option->speed = 11;
}
