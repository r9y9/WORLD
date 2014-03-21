//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise. All Rights Reserved.
// Author: morise [at] fc.ritsumei.ac.jp (Masanori Morise)
//
// Spectral envelope estimation based on STAR (Synchronous Technique and Adroit
// Restoration).
// Please see styleguide.txt to show special rules on names of variables
// and fnctions.
//-----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include "./star.h"
#include "./matlabfunctions.h"
#include "./constant_numbers.h"

namespace {

//-----------------------------------------------------------------------------
// AdroitSmoothing() carries out the spectral smoothing by rectangular window
// whose length is F0.
// This function is only used in StarGeneralBody().
//-----------------------------------------------------------------------------
void AdroitSmoothing(double current_f0, int fs, int fft_size,
    double *power_spectrum, double *star_spectrum) {
  int boundary = static_cast<int>(current_f0 /
    (static_cast<double>(fs) / fft_size)) + 1;
  double *mirroring_spectrum = new double[fft_size + boundary * 2 + 1];

  for (int i = 0; i < boundary; ++i)
    mirroring_spectrum[i] = power_spectrum[boundary - i];
  for (int i = boundary; i < fft_size / 2 + boundary; ++i)
    mirroring_spectrum[i] = power_spectrum[i - boundary];
  for (int i = fft_size / 2 + boundary;
      i < fft_size / 2 + boundary * 2 + 1; ++i)
    mirroring_spectrum[i] =
    power_spectrum[fft_size / 2 - (i - (fft_size / 2 + boundary)) - 1];

  int tmp = static_cast<int>(current_f0 * fft_size / fs);

  double *mirroring_segment = new double[fft_size * 2];
  mirroring_segment[0] = log(mirroring_spectrum[0]) * fs / fft_size;
  for (int i = 1; i < fft_size / 2 + boundary * 2 + 1; ++i)
    mirroring_segment[i] = log(mirroring_spectrum[i]) * fs / fft_size +
    mirroring_segment[i - 1];

  double *frequency_axis = new double[fft_size / 2 + 1];
  for (int i = 0; i <= fft_size / 2; ++i)
    frequency_axis[i] = static_cast<double>(i) / fft_size *
    fs - current_f0 / 2.0;

  double *low_levels = new double[fft_size / 2 + 1];
  double *high_levels = new double[fft_size / 2 + 1];
  double origin_of_mirroring_axis =
    -(static_cast<double>(boundary) - 0.5) * fs / fft_size;
  double discrete_frequency_interval = static_cast<double>(fs) / fft_size;

  interp1Q(origin_of_mirroring_axis, discrete_frequency_interval,
      mirroring_segment, fft_size / 2 + boundary * 2 + 1, frequency_axis,
      fft_size / 2 + 1, low_levels);
  for (int i = 0; i <= fft_size / 2; ++i) frequency_axis[i] += current_f0;

  interp1Q(origin_of_mirroring_axis, discrete_frequency_interval,
      mirroring_segment, fft_size / 2 + boundary * 2 + 1, frequency_axis,
      fft_size / 2 + 1, high_levels);

  for (int i = 0; i <= fft_size / 2; ++i)
    star_spectrum[i] = exp((high_levels[i] - low_levels[i]) / current_f0);

  delete[] low_levels;
  delete[] high_levels;
  delete[] mirroring_segment;
  delete[] frequency_axis;
  delete[] mirroring_spectrum;
}

//-----------------------------------------------------------------------------
// GetPowerSpectrum() carries out (1) designing the window,
// (2) windowing the waveform and (3) calculation of the power_spectrum
//-----------------------------------------------------------------------------
void GetPowerSpectrum(double *x, int x_length, int fs, double current_f0,
    double temporal_position, ForwardRealFFT *forward_real_fft,
    double *power_spectrum) {
  int half_window_length = matlab_round(3.0 * fs / current_f0 / 2.0);
  int *base_index = new int[half_window_length * 2 + 1];
  int *index = new int[half_window_length * 2 + 1];

  for (int i = -half_window_length; i <= half_window_length; ++i)
    base_index[i + half_window_length] = i;
  for (int i = 0; i <= half_window_length * 2; ++i)
    index[i] = std::min(x_length, std::max(1,
    matlab_round(temporal_position * fs + 1 + base_index[i]))) - 1;

  // Designing of the window function
  double *window  = new double[half_window_length * 2 + 1];
  double average = 0.0;
  double position;
  for (int i = 0; i <= half_window_length * 2; ++i) {
    position = static_cast<double>(base_index[i]) / fs / (3.0 / 2.0) +
      (temporal_position * fs - matlab_round(temporal_position * fs)) / fs;
    window[i] = 0.5 * cos(world::kPi * position * current_f0) + 0.5;
    average += window[i] * window[i];
  }
  average = sqrt(average);
  for (int i = 0; i <= half_window_length * 2; ++i)
    window[i] /= average;

  // Windowing and FFT
  for (int i = 0; i <= half_window_length * 2; ++i)
    forward_real_fft->waveform[i] = x[index[i]] * window[i];
  for (int i = half_window_length * 2 + 1; i < forward_real_fft->fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;
  fft_execute(forward_real_fft->forward_fft);

  // Calculation of the power spectrum.
  for (int i = 1; i <= forward_real_fft->fft_size / 2; ++i)
    power_spectrum[i] = forward_real_fft->spectrum[i][0] *
    forward_real_fft->spectrum[i][0] +
    forward_real_fft->spectrum[i][1] *
    forward_real_fft->spectrum[i][1];
  power_spectrum[0] = power_spectrum[1];

  delete[] window;
  delete[] base_index;
  delete[] index;
}

//-----------------------------------------------------------------------------
// StarGeneralBody() calculates a spectral envelope at a temporal position.
// This function is only used in Star().
// Caution:
//   windowed_waveform, y_spectrum and forward_fft is allocated in advance in
//   Star() to speed up the processing. If you want to develop real-time
//   application, you should modify this function not to use these arguments
//   and edit this function.
//-----------------------------------------------------------------------------
void StarGeneralBody(double *x, int x_length, int fs, double current_f0,
    double temporal_position, ForwardRealFFT *forward_real_fft,
    double * star_spectrum) {
  double *power_spectrum = new double[forward_real_fft->fft_size];

  // Synchronous analysis
  GetPowerSpectrum(x, x_length, fs, current_f0, temporal_position,
      forward_real_fft, power_spectrum);

  // Adroit smoothing
  AdroitSmoothing(current_f0, fs, forward_real_fft->fft_size,
      power_spectrum, star_spectrum);

  delete[] power_spectrum;
}

}  // namespace

void Star(double *x, int x_length, int fs, double *time_axis, double *f0,
    int f0_length, double **spectrogram) {
  double frame_period = (time_axis[1] - time_axis[0]) * 1000.0;

  int fft_size = GetFFTSizeForStar(fs);

  double *star_spectrum = new double[fft_size];

  // Following three variables are shared in StarGeneralBody()
  ForwardRealFFT forward_real_fft = {0};
  InitializeForwardRealFFT(fft_size, &forward_real_fft);

  double current_f0;
  for (int i = 0; i < f0_length; ++i) {
    current_f0 = f0[i] <= world::kFloorF0 ? world::kDefaultF0 : f0[i];
    StarGeneralBody(x, x_length, fs, current_f0, time_axis[i],
        &forward_real_fft, star_spectrum);
    for (int j = 0; j <= fft_size / 2; ++j)
      spectrogram[i][j] = star_spectrum[j];
  }

  DestroyForwardRealFFT(&forward_real_fft);
  delete[] star_spectrum;
}

int GetFFTSizeForStar(int fs) {
  return static_cast<int>(pow(2.0, 1.0 +
    static_cast<int>(log(3.0 * fs / world::kFloorF0 + 1) / world::kLog2)));
}
