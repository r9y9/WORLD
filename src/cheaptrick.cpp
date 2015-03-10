//-----------------------------------------------------------------------------
// Copyright 2012-2014 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//
// Spectral envelope estimation on the basis of the idea of CheapTrick.
//-----------------------------------------------------------------------------
#include "./cheaptrick.h"

#include <math.h>

#include "./constantnumbers.h"
#include "./matlabfunctions.h"

namespace {

//-----------------------------------------------------------------------------
// SmoothingWithRecovery() carries out the spectral smoothing and spectral
// recovery on the Cepstrum domain.
//-----------------------------------------------------------------------------
void SmoothingWithRecovery(double current_f0, int fs, int fft_size,
    ForwardRealFFT *forward_real_fft, InverseRealFFT *inverse_real_fft,
    double *spectral_envelope) {
  const double q1 = -0.09;  // Please see the reference in CheapTrick.
  double *smoothing_lifter = new double[fft_size];
  double *compensation_lifter = new double[fft_size];

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

  for (int i = 1; i < fft_size / 2; ++i)
    forward_real_fft->waveform[fft_size - i] = forward_real_fft->waveform[i];
  fft_execute(forward_real_fft->forward_fft);

  for (int i = 0; i <= fft_size / 2; ++i) {
    inverse_real_fft->spectrum[i][0] = forward_real_fft->spectrum[i][0] *
      smoothing_lifter[i] * compensation_lifter[i] / fft_size;
    inverse_real_fft->spectrum[i][1] = 0.0;
  }
  fft_execute(inverse_real_fft->inverse_fft);

  for (int i = 0; i <= fft_size / 2; ++i)
    spectral_envelope[i] = exp(inverse_real_fft->waveform[i]);

  delete[] smoothing_lifter;
  delete[] compensation_lifter;
}

//-----------------------------------------------------------------------------
// SetParametersForLinearSmoothing()
//-----------------------------------------------------------------------------
void SetParametersForLinearSmoothing(int boundary, int fft_size, int fs,
    double current_f0, double *power_spectrum, double *mirroring_spectrum,
    double *mirroring_segment, double *frequency_axis) {
  for (int i = 0; i < boundary; ++i)
    mirroring_spectrum[i] = power_spectrum[boundary - i];
  for (int i = boundary; i < fft_size / 2 + boundary; ++i)
    mirroring_spectrum[i] = power_spectrum[i - boundary];
  for (int i = fft_size / 2 + boundary; i <= fft_size / 2 + boundary * 2; ++i)
    mirroring_spectrum[i] =
      power_spectrum[fft_size / 2 - (i - (fft_size / 2 + boundary))];

  mirroring_segment[0] = mirroring_spectrum[0] * fs / fft_size;
  for (int i = 1; i < fft_size / 2 + boundary * 2 + 1; ++i)
    mirroring_segment[i] = mirroring_spectrum[i] * fs / fft_size +
    mirroring_segment[i - 1];

  for (int i = 0; i <= fft_size / 2; ++i)
    frequency_axis[i] = static_cast<double>(i) / fft_size *
    fs - current_f0 / 3.0;
}

//-----------------------------------------------------------------------------
// LinearSmoothing() carries out the spectral smoothing by rectangular window
// whose length is 2F0 / 3.
// Note that the output is "Logarithmic" power spectrum.
//-----------------------------------------------------------------------------
void LinearSmoothing(double current_f0, int fs, int fft_size,
    ForwardRealFFT *forward_real_fft) {
  int boundary = static_cast<int>(current_f0 * fft_size / fs) + 1;

  // These parameters are set by the other function.
  double *power_spectrum = forward_real_fft->waveform;
  double *mirroring_spectrum = new double[fft_size / 2 + boundary * 2 + 1];
  double *mirroring_segment = new double[fft_size / 2 + boundary * 2 + 1];
  double *frequency_axis = new double[fft_size / 2 + 1];
  SetParametersForLinearSmoothing(boundary, fft_size, fs, current_f0,
      power_spectrum, mirroring_spectrum, mirroring_segment, frequency_axis);

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

  double *smoothed_spectrum = forward_real_fft->waveform;
  for (int i = 0; i <= fft_size / 2; ++i)
    smoothed_spectrum[i] =
      log((high_levels[i] - low_levels[i]) * 1.5 / current_f0);

  delete[] mirroring_spectrum;
  delete[] mirroring_segment;
  delete[] frequency_axis;
  delete[] low_levels;
  delete[] high_levels;
}

//-----------------------------------------------------------------------------
// GetPowerSpectrum() calculates the power_spectrum with DC correction
//-----------------------------------------------------------------------------
void GetPowerSpectrum(int fs, double current_f0, int fft_size,
    ForwardRealFFT *forward_real_fft) {
  int half_window_length = matlab_round(1.5 * fs / current_f0);

  // FFT
  for (int i = half_window_length * 2 + 1; i < fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;
  fft_execute(forward_real_fft->forward_fft);

  // Calculation of the power spectrum.
  double *power_spectrum = forward_real_fft->waveform;
  for (int i = 0; i <= fft_size / 2; ++i)
    power_spectrum[i] =
      forward_real_fft->spectrum[i][0] * forward_real_fft->spectrum[i][0] +
      forward_real_fft->spectrum[i][1] * forward_real_fft->spectrum[i][1];

  // DC correction
  int upper_limit = 1 +
    static_cast<int>(1.2 * current_f0 * fft_size / fs);
  double *low_frequency_replica = new double[upper_limit];
  double *low_frequency_axis = new double[upper_limit];

  for (int i = 0; i < upper_limit; ++i)
    low_frequency_axis[i] = static_cast<double>(i) * fs / fft_size;

  // Bug fix! 2014/10/11 by M. Morise
  int upper_limit_replica = 1 + static_cast<int>(current_f0 * fft_size / fs);
  interp1Q(current_f0 - low_frequency_axis[0],
      -static_cast<double>(fs) / fft_size, power_spectrum, upper_limit + 1,
      low_frequency_axis, upper_limit_replica, low_frequency_replica);

  for (int i = 0; i < upper_limit_replica; ++i)
    power_spectrum[i] += low_frequency_replica[i];

  delete[] low_frequency_replica;
  delete[] low_frequency_axis;
}

//-----------------------------------------------------------------------------
// SetparametersForGetWindowedWaveform()
//-----------------------------------------------------------------------------
void SetparametersForGetWindowedWaveform(int half_window_length, int x_length,
    double temporal_position, int fs, double current_f0, int *base_index,
    int *index, double *window) {
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
}

//-----------------------------------------------------------------------------
// GetWindowedWaveform() windows the waveform by pitch synchronous window
//-----------------------------------------------------------------------------
void GetWindowedWaveform(double *x, int x_length, int fs, double current_f0,
    double temporal_position, ForwardRealFFT *forward_real_fft) {
  int half_window_length = matlab_round(1.5 * fs / current_f0);

  int *base_index = new int[half_window_length * 2 + 1];
  int *index = new int[half_window_length * 2 + 1];
  double *window  = new double[half_window_length * 2 + 1];

  SetparametersForGetWindowedWaveform(half_window_length, x_length,
      temporal_position, fs, current_f0, base_index, index, window);

  // Pitch synchronous windowing
  double *waveform = forward_real_fft->waveform;
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
//   windowed_waveform, y_spectrum and forward_fft is allocated in advance
//   to speed up the processing. If you want to develop real-time
//   application, you should modify this function not to use these arguments
//   and edit this function.
//-----------------------------------------------------------------------------
void CheapTrickGeneralBody(double *x, int x_length, int fs, double current_f0,
    int fft_size, double temporal_position, ForwardRealFFT *forward_real_fft,
    InverseRealFFT *inverse_real_fft, double *spectral_envelope) {
  // Synchronous windowing
  GetWindowedWaveform(x, x_length, fs, current_f0, temporal_position,
      forward_real_fft);

  // Calculate power spectrum with DC correction
  // Note: The calculated power spectrum is stored in an array for waveform.
  // In this imprementation, power spectrum is transformed by FFT (NOT IFFT).
  // However, the same result is obtained.
  // This is tricky but important for simple implementation.
  GetPowerSpectrum(fs, current_f0, fft_size, forward_real_fft);

  // Smoothing of the power (linear axis)
  // "work_space->forward_real_fft.waveform" is the power spectrum.
  LinearSmoothing(current_f0, fs, fft_size, forward_real_fft);

  // Smoothing and spectral recovery on the cepstrum domain.
  SmoothingWithRecovery(current_f0, fs, fft_size, forward_real_fft,
      inverse_real_fft, spectral_envelope);
}

}  // namespace


DLLEXPORT int GetFFTSizeForCheapTrick(int fs) {
  return static_cast<int>(pow(2.0, 1.0 +
      static_cast<int>(log(3.0 * fs / world::kFloorF0 + 1) / world::kLog2)));
}

DLLEXPORT void CheapTrick(double *x, int x_length, int fs, double *time_axis,
    double *f0, int f0_length, double **spectrogram) {
  int fft_size = GetFFTSizeForCheapTrick(fs);
  double *spectral_envelope = new double[fft_size];

  ForwardRealFFT forward_real_fft = {0};
  InitializeForwardRealFFT(fft_size, &forward_real_fft);
  InverseRealFFT inverse_real_fft = {0};
  InitializeInverseRealFFT(fft_size, &inverse_real_fft);

  double current_f0;
  for (int i = 0; i < f0_length; ++i) {
    current_f0 = f0[i] <= world::kFloorF0 ? world::kDefaultF0 : f0[i];
    CheapTrickGeneralBody(x, x_length, fs, current_f0, fft_size,
        time_axis[i], &forward_real_fft, &inverse_real_fft, spectral_envelope);
    for (int j = 0; j <= fft_size / 2; ++j)
      spectrogram[i][j] = spectral_envelope[j];
  }

  DestroyForwardRealFFT(&forward_real_fft);
  DestroyInverseRealFFT(&inverse_real_fft);
  delete[] spectral_envelope;
}
