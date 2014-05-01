//-----------------------------------------------------------------------------
// Copyright 2012-2014 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//
// Spectral envelope estimation based on CheapTrick.
//-----------------------------------------------------------------------------
#include "./cheaptrick.h"

#include <math.h>
#include <stdlib.h>

#include "./constantnumbers.h"
#include "./matlabfunctions.h"

/*
メモ：f0に依存せず長さが固定の配列は，毎回確保，解放するのではなく
最初に作業用スペースとして確保し，それを使いまわすほうが効果的．
forward_real_fft
inverse_real_fft

smoothing_lifter
compensation_lifter

frequency_axis
low_level
high_level

power_spectrum
smoothed_spectrum
waveform

*/

namespace {

//-----------------------------------------------------------------------------
// SmoothingWithRecovery() carries out the spectral smoothing and spectral
// recovery on the Cepstrum domain.
//-----------------------------------------------------------------------------
void SmoothingWithRecovery(double current_f0, int fs,
    ForwardRealFFT *forward_real_fft, InverseRealFFT *inverse_real_fft,
    double *smoothed_spectrum, double *spectral_envelope) {
  const double q1 = -0.09; // Please see the reference in CheapTrick.
  double *smoothing_lifter = new double[forward_real_fft->fft_size];
  double *compensation_lifter = new double[forward_real_fft->fft_size];

  smoothing_lifter[0] = 1;
  compensation_lifter[0] = (1.0 - 2.0 * q1) + 2.0 * q1;
  double quefrency;
  for (int i = 1; i <= forward_real_fft->fft_size / 2; ++i) {
    quefrency = static_cast<double>(i) / fs;
    smoothing_lifter[i] = sin(world::kPi * current_f0 * quefrency) /
      (world::kPi * current_f0 * quefrency);
    compensation_lifter[i] = (1.0 - 2.0 * q1) + 2.0 * q1 *
      cos(2.0 * world::kPi * quefrency * current_f0);
  }

  forward_real_fft->waveform[0] = log(smoothed_spectrum[0]);
  forward_real_fft->waveform[forward_real_fft->fft_size / 2] =
    log(smoothed_spectrum[forward_real_fft->fft_size / 2]);
  for (int i = 1; i < forward_real_fft->fft_size / 2; ++i) {
    forward_real_fft->waveform[i] = log(smoothed_spectrum[i]);
    forward_real_fft->waveform[forward_real_fft->fft_size - i] =
      forward_real_fft->waveform[i];
  }
  fft_execute(forward_real_fft->forward_fft);

  for (int i = 0; i <= inverse_real_fft->fft_size / 2; ++i) {
    inverse_real_fft->spectrum[i][0] = forward_real_fft->spectrum[i][0] *
      smoothing_lifter[i] * compensation_lifter[i] / inverse_real_fft->fft_size;
    inverse_real_fft->spectrum[i][1] = 0.0;
  }
  fft_execute(inverse_real_fft->inverse_fft);

  for (int i = 0; i <= inverse_real_fft->fft_size / 2; ++i)
    spectral_envelope[i] = exp(inverse_real_fft->waveform[i]);

  delete[] smoothing_lifter;
  delete[] compensation_lifter;
}

//-----------------------------------------------------------------------------
// LinearSmoothing() carries out the spectral smoothing by rectangular window
// whose length is 2F0 / 3.
//-----------------------------------------------------------------------------
void LinearSmoothing(double current_f0, int fs, int fft_size,
    double *power_spectrum, double *smoothed_spectrum) {
  int boundary = static_cast<int>(current_f0 /
    (static_cast<double>(fs) / fft_size)) + 1;
  double *mirroring_spectrum = new double[fft_size / 2 + boundary * 2 + 1];

  for (int i = 0; i < boundary; ++i)
    mirroring_spectrum[i] = power_spectrum[boundary - i];
  for (int i = boundary; i < fft_size / 2 + boundary; ++i)
    mirroring_spectrum[i] = power_spectrum[i - boundary];
  for (int i = fft_size / 2 + boundary;
      i < fft_size / 2 + boundary * 2 + 1; ++i)
    mirroring_spectrum[i] =
    power_spectrum[fft_size / 2 - (i - (fft_size / 2 + boundary))];

  double *mirroring_segment = new double[fft_size / 2+ boundary * 2 + 1];
  mirroring_segment[0] = mirroring_spectrum[0] * fs / fft_size;
  for (int i = 1; i < fft_size / 2 + boundary * 2 + 1; ++i)
    mirroring_segment[i] = mirroring_spectrum[i] * fs / fft_size +
    mirroring_segment[i - 1];

  double *frequency_axis = new double[fft_size / 2 + 1];
  for (int i = 0; i <= fft_size / 2; ++i)
    frequency_axis[i] = static_cast<double>(i) / fft_size *
    fs - current_f0 / 3.0;

  double *low_levels = new double[fft_size / 2 + 1];
  double *high_levels = new double[fft_size / 2 + 1];
  double origin_of_mirroring_axis =
    -(static_cast<double>(boundary) - 0.5) * fs / fft_size;
  double discrete_frequency_interval = static_cast<double>(fs) / fft_size;

  interp1Q(origin_of_mirroring_axis, discrete_frequency_interval,
      mirroring_segment, fft_size / 2 + boundary * 2 + 1, frequency_axis,
      fft_size / 2 + 1, low_levels);
  for (int i = 0; i <= fft_size / 2; ++i)
    frequency_axis[i] += 2.0 * current_f0 / 3.0;

  interp1Q(origin_of_mirroring_axis, discrete_frequency_interval,
      mirroring_segment, fft_size / 2 + boundary * 2 + 1, frequency_axis,
      fft_size / 2 + 1, high_levels);

  for (int i = 0; i <= fft_size / 2; ++i)
    smoothed_spectrum[i] = (high_levels[i] - low_levels[i]) * 1.5 / current_f0;

  delete[] mirroring_spectrum;
  delete[] mirroring_segment;
  delete[] frequency_axis;
  delete[] low_levels;
  delete[] high_levels;
}



//-----------------------------------------------------------------------------
// GetPowerSpectrum() calculates the power_spectrum with DC correction
//-----------------------------------------------------------------------------
void GetPowerSpectrum(double *waveform, int fs, double current_f0,
    ForwardRealFFT *forward_real_fft, double *power_spectrum) {
  int half_window_length = matlab_round(3.0 * fs / current_f0 / 2.0);

  // Windowing and FFT
  for (int i = 0; i <= half_window_length * 2; ++i)
    forward_real_fft->waveform[i] = waveform[i];
  for (int i = half_window_length * 2 + 1; i < forward_real_fft->fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;
  fft_execute(forward_real_fft->forward_fft);

  // Calculation of the power spectrum.
  for (int i = 0; i <= forward_real_fft->fft_size / 2; ++i)
    power_spectrum[i] = forward_real_fft->spectrum[i][0] *
      forward_real_fft->spectrum[i][0] +
      forward_real_fft->spectrum[i][1] *
      forward_real_fft->spectrum[i][1];

  // DC correction
  int upper_limit = 1 + 
    static_cast<int>(1.2 * current_f0 / (fs / forward_real_fft->fft_size));
  double *low_frequency_replica = new double[upper_limit];
  double *low_frequency_axis = new double[upper_limit];
  double *reverse_frequency_axis = new double[upper_limit + 1];
  double *dammy_power_spectrum = new double[upper_limit + 1];

  for (int i = 0; i < upper_limit; ++i) {
    dammy_power_spectrum[i] = power_spectrum[upper_limit - i - 1];
    low_frequency_axis[i] = static_cast<double>(i) *
      fs / forward_real_fft->fft_size;
    reverse_frequency_axis[upper_limit - i - 1] =
      current_f0 - low_frequency_axis[i];
  }
  reverse_frequency_axis[upper_limit] =
    reverse_frequency_axis[upper_limit - 1] * 2 -
    reverse_frequency_axis[upper_limit - 2];
  dammy_power_spectrum[upper_limit] =
    dammy_power_spectrum[upper_limit - 1] * 2 -
    dammy_power_spectrum[upper_limit - 2];

  interp1(reverse_frequency_axis, dammy_power_spectrum, upper_limit + 1,
    low_frequency_axis, upper_limit, low_frequency_replica);

  upper_limit = 1 + 
    static_cast<int>(current_f0 / (fs / forward_real_fft->fft_size));
  for (int i = 0; i < upper_limit; ++i)
    power_spectrum[i] += low_frequency_replica[i];

  delete[] low_frequency_replica;
  delete[] low_frequency_axis;
  delete[] reverse_frequency_axis;
  delete[] dammy_power_spectrum;
}

//-----------------------------------------------------------------------------
// GetWindowedWaveform() windows the waveform by pitch synchronous window
//-----------------------------------------------------------------------------
void GetWindowedWaveform(double *x, int x_length, int fs, double current_f0,
    double temporal_position, double *waveform) {
  int half_window_length = matlab_round(1.5 * fs / current_f0);
  int *base_index = new int[half_window_length * 2 + 1];
  int *index = new int[half_window_length * 2 + 1];
  double *window  = new double[half_window_length * 2 + 1];

  for (int i = -half_window_length; i <= half_window_length; ++i)
    base_index[i + half_window_length] = i;
  for (int i = 0; i <= half_window_length * 2; ++i)
    index[i] = MyMin(x_length - 1, MyMax(0,
        matlab_round(temporal_position * fs + base_index[i])));

  // Designing of the window function
  double average = 0.0;
  double position;
  double bias = temporal_position * fs - matlab_round(temporal_position * fs);
  for (int i = 0; i <= half_window_length * 2; ++i) {
    position = (static_cast<double>(base_index[i]) / 1.5 + bias) / fs;
    window[i] = 0.5 * cos(world::kPi * position * current_f0) + 0.5;
    average += window[i] * window[i];
  }
  average = sqrt(average);
  for (int i = 0; i <= half_window_length * 2; ++i) window[i] /= average;

  // Pitch synchronous windowing
  for (int i = 0; i <= half_window_length * 2; ++i)
    waveform[i] = x[index[i]] * window[i] + randn() * 0.000000000000001;
  double tmp_weight1 = 0;
  double tmp_weight2 = 0;
  for (int i = 0; i <= half_window_length * 2; ++i) {
    tmp_weight1 += waveform[i];
    tmp_weight2 += window[i];
  }
  double weighting_coefficient = tmp_weight1 / tmp_weight2;
  for (int i = 0; i <= half_window_length * 2; ++i)
    waveform[i] -= window[i] * weighting_coefficient; 

  delete[] base_index;
  delete[] index;
  delete[] window;
}

//-----------------------------------------------------------------------------
// CheapTrickGeneralBody() calculates a spectral envelope at a temporal 
// position. This function is only used in CheapTrick().
// Caution:
//   windowed_waveform, y_spectrum and forward_fft is allocated in advance in
//   Star() to speed up the processing. If you want to develop real-time
//   application, you should modify this function not to use these arguments
//   and edit this function.
//-----------------------------------------------------------------------------
void CheapTrickGeneralBody(double *x, int x_length, int fs, double current_f0,
    double temporal_position, ForwardRealFFT *forward_real_fft,
    InverseRealFFT *inverse_real_fft, double * spectral_envelope) {
  double *power_spectrum = new double[forward_real_fft->fft_size];
  double *smoothed_spectrum = new double[forward_real_fft->fft_size];
  double *waveform = new double[forward_real_fft->fft_size];

  // Synchronous windowing
  GetWindowedWaveform(x, x_length, fs, current_f0, temporal_position,
      waveform);

  // Calculate power spectrum with DC correction
  GetPowerSpectrum(waveform, fs, current_f0, forward_real_fft, power_spectrum);

  // Smoothing of the power (linear axis)
  LinearSmoothing(current_f0, fs, forward_real_fft->fft_size,
      power_spectrum, smoothed_spectrum);

  SmoothingWithRecovery(current_f0, fs, forward_real_fft, inverse_real_fft,
    smoothed_spectrum, spectral_envelope);

  delete[] power_spectrum;
  delete[] smoothed_spectrum;
  delete[] waveform;
}

}  // namespace


int GetFFTSizeForCheapTrick(int fs) {
  return static_cast<int>(pow(2.0, 1.0 +
      static_cast<int>(log(3.0 * fs / world::kFloorF0 + 1) / world::kLog2)));
}

void CheapTrick(double *x, int x_length, int fs, double *time_axis, double *f0,
    int f0_length, double **spectrogram) {
  int fft_size = GetFFTSizeForCheapTrick(fs);

  double *spectral_envelope = new double[fft_size];

  // Following three variables are shared in StarGeneralBody()
  ForwardRealFFT forward_real_fft = {0};
  InitializeForwardRealFFT(fft_size, &forward_real_fft);
  InverseRealFFT inverse_real_fft = {0};
  InitializeInverseRealFFT(fft_size, &inverse_real_fft);

  double current_f0;
  for (int i = 0; i < f0_length; ++i) {
    current_f0 = f0[i] <= world::kFloorF0 ? world::kDefaultF0 : f0[i];
    CheapTrickGeneralBody(x, x_length, fs, current_f0, time_axis[i],
        &forward_real_fft, &inverse_real_fft, spectral_envelope);
    for (int j = 0; j <= fft_size / 2; ++j)
      spectrogram[i][j] = spectral_envelope[j];
  }

  DestroyForwardRealFFT(&forward_real_fft);
  DestroyInverseRealFFT(&inverse_real_fft);
  delete[] spectral_envelope;
}
