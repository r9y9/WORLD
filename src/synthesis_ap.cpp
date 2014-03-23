//-----------------------------------------------------------------------------
// Copyright 2012-2013 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//
// Voice synthesis based on f0, spectrogram and aperiodicity.
// forward_real_fft, inverse_real_fft and minimum_phase are used to speed up.
//-----------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "./common.h"
#include "./constantnumbers.h"
#include "./matlabfunctions.h"
#include "./synthesis_ap.h"
#include "./tandem_ap.h"

namespace {

//-----------------------------------------------------------------------------
// GetGlottalPulse() calculates the glottal pulse based on periodic response
// and aperiodic response.
//-----------------------------------------------------------------------------
void GetGlottalPulse(double f0, int fft_size, double *periodic_response,
    double *aperiodic_response, int noise_size, double *y) {
  if (f0 != 0) {
    for (int i = 0; i < fft_size; ++i)
      y[i] = periodic_response[i] * sqrt(static_cast<double>(noise_size)) +
      aperiodic_response[i];
  } else {
    for (int i = 0; i < fft_size; ++i)
      y[i] = aperiodic_response[i];
  }
  for (int i = 0; i < fft_size; ++i) y[i] /= fft_size;
}

//-----------------------------------------------------------------------------
// CalculateAperiodicity() transforms the input aperiodicity in each band
// into the aperiodicity spectrum whose length is fft_size.
//-----------------------------------------------------------------------------
void CalculateAperiodicity(double *aperiodicity, int number_of_bands,
    int fft_size, double f0, int fs, double target_f0, double *periodic_spec) {
  if (f0 == 0) {
    for (int i = 0; i <= fft_size / 2; ++i) periodic_spec[i] = 0.0;
    return;
  }
  double *ap = new double[number_of_bands + 1];
  double *axis = new double[number_of_bands + 1];
  double *w = new double[fft_size / 2 + 1];
  double *tmp_ap = new double[fft_size / 2 + 1];

  double *cutoff_list = new double[number_of_bands];
  for (int i = 0; i < number_of_bands; ++i)
    cutoff_list[i] = fs / pow(2.0, i + 2.0);

  const double kMySafeGuardLogMinimum = -27.631021115928547;
  ap[0] = kMySafeGuardLogMinimum;
  axis[0] = 0.0;
  for (int i = 0; i < number_of_bands - 1; ++i) {
    ap[i + 1] = log(aperiodicity[i]);
    axis[i + 1] = cutoff_list[number_of_bands - i - 2];
  }
  ap[number_of_bands] = log(aperiodicity[number_of_bands - 1]);
  axis[number_of_bands] = fs / 2.0;

  double stretching_factor = MyMax(f0, target_f0) / target_f0;
  for (int i = 0; i <= fft_size / 2; ++i)
    w[i] = static_cast<double>(i * fs) / fft_size;
  interp1(axis, ap, number_of_bands + 1, w, fft_size / 2 + 1, tmp_ap);
  for (int i = 0; i < number_of_bands - 1; ++i)
    axis[i + 1] *= stretching_factor;
  axis[number_of_bands] = fs / 2.0 * stretching_factor;

  interp1(axis, ap, number_of_bands + 1, w, fft_size / 2 + 1, periodic_spec);

  for (int i = 0; i <= fft_size / 2; ++i)
    periodic_spec[i] = 1.0 -
    MyMin(exp(tmp_ap[i] * 2.0), exp(periodic_spec[i] * 2.0));

  delete[] tmp_ap;
  delete[] cutoff_list;
  delete[] w;
  delete[] axis;
  delete[] ap;
}

//-----------------------------------------------------------------------------
// GetAperiodicResponse() calculates an aperiodic response.
//-----------------------------------------------------------------------------
void GetAperiodicResponse(int noise_size, int fft_size,
    double *spectrum, double *aperiodic_ratio,
    ForwardRealFFT *forward_real_fft, InverseRealFFT *inverse_real_fft,
    MinimumPhaseAnalysis *minimum_phase,
    double *aperiodic_response) {
  for (int i = 0; i < noise_size; ++i) forward_real_fft->waveform[i] = randn();
  for (int i = noise_size; i < fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;
  fft_execute(forward_real_fft->forward_fft);

  for (int i = 0; i <= minimum_phase->fft_size / 2; ++i)
    minimum_phase->log_spectrum[i] =
    log(spectrum[i] *
    ((1 - aperiodic_ratio[i]) + world::kMySafeGuardMinimum)) / 2.0;
  GetMinimumPhaseSpectrum(minimum_phase);

  for (int i = 0; i <= fft_size / 2; ++i) {
    inverse_real_fft->spectrum[i][0] =
      minimum_phase->minimum_phase_spectrum[i][0] *
      forward_real_fft->spectrum[i][0] -
      minimum_phase->minimum_phase_spectrum[i][1] *
      forward_real_fft->spectrum[i][1];
    inverse_real_fft->spectrum[i][1] =
      minimum_phase->minimum_phase_spectrum[i][0] *
      forward_real_fft->spectrum[i][1] +
      minimum_phase->minimum_phase_spectrum[i][1] *
      forward_real_fft->spectrum[i][0];
  }
  fft_execute(inverse_real_fft->inverse_fft);
  for (int i = 0; i < fft_size; ++i)
    aperiodic_response[i] = inverse_real_fft->waveform[i];
}

//-----------------------------------------------------------------------------
// GetPeriodicResponse() calculates an aperiodic response.
//-----------------------------------------------------------------------------
void GetPeriodicResponse(int fft_size, double *spectrum,
    double *aperiodic_ratio, InverseRealFFT *inverse_real_fft,
    MinimumPhaseAnalysis *minimum_phase,
    double *periodic_response) {
  for (int i = 0; i <= minimum_phase->fft_size / 2; ++i)
    minimum_phase->log_spectrum[i] =
    log(spectrum[i] * aperiodic_ratio[i]) / 2.0;
  GetMinimumPhaseSpectrum(minimum_phase);
  for (int i = 0; i <= fft_size / 2; ++i) {
    inverse_real_fft->spectrum[i][0] =
      minimum_phase->minimum_phase_spectrum[i][0];
    inverse_real_fft->spectrum[i][1] =
      minimum_phase->minimum_phase_spectrum[i][1];
  }
  fft_execute(inverse_real_fft->inverse_fft);
  for (int i = 0; i < fft_size; ++i)
    periodic_response[i] = inverse_real_fft->waveform[i];
}

//-----------------------------------------------------------------------------
// GetOneFrameSegment() calculates a glottal pulse at a time.
//-----------------------------------------------------------------------------
void GetOneFrameSegment(double *f0, double **spectrogram, int fft_size,
    double **aperiodicity, int number_of_bands, double target_f0,
    double frame_period, double current_time, int fs, double default_f0,
    ForwardRealFFT *forward_real_fft, InverseRealFFT *inverse_real_fft,
    MinimumPhaseAnalysis *minimum_phase, double *y) {
  double *aperiodic_ratio = new double[fft_size];
  double *aperiodic_response = new double[fft_size];
  double *periodic_response = new double[fft_size];

  int current_frame = matlab_round(current_time / (frame_period / 1000.0));
  int noise_size = static_cast<int>((current_time + 1.0 /
    (f0[current_frame] == 0.0 ? default_f0 : f0[current_frame])) * fs) -
    static_cast<int>(current_time * fs);

  // Calculation of the aperiodicity at each discrete frequency
  CalculateAperiodicity(aperiodicity[current_frame], number_of_bands,
      fft_size, f0[current_frame], fs, target_f0, aperiodic_ratio);

  // Synthesis of the aperiodic response
  GetAperiodicResponse(noise_size, fft_size, spectrogram[current_frame],
      aperiodic_ratio, forward_real_fft, inverse_real_fft, minimum_phase,
      aperiodic_response);

  // Synthesis of the periodic response.
  // If f0 is zero, we cannot synthesize it.
  if (f0[current_frame] != 0) {
    GetPeriodicResponse(fft_size, spectrogram[current_frame],
        aperiodic_ratio, inverse_real_fft, minimum_phase,
        periodic_response);
  }

  GetGlottalPulse(f0[current_frame], fft_size, periodic_response,
      aperiodic_response, noise_size, y);

  delete[] periodic_response;
  delete[] aperiodic_response;
  delete[] aperiodic_ratio;
}

}  // namespace

void SynthesisFromAperiodicityOld(double *f0, int f0_length, double **spectrogram,
    int fft_size, double **aperiodicity, int number_of_bands, double target_f0,
    double frame_period, int fs, int y_length, double *y) {
  double *impulse_response = new double[fft_size];

  for (int i = 0; i < y_length; ++i) y[i] = 0.0;

  MinimumPhaseAnalysis minimum_phase = {0};
  InitializeMinimumPhaseAnalysis(fft_size, &minimum_phase);
  InverseRealFFT inverse_real_fft = {0};
  InitializeInverseRealFFT(fft_size, &inverse_real_fft);
  ForwardRealFFT forward_real_fft = {0};
  InitializeForwardRealFFT(fft_size, &forward_real_fft);

  double current_time = 0.0;
  int current_position = 0;
  int current_frame = 0;
  while (1) {
    GetOneFrameSegment(f0, spectrogram, fft_size, aperiodicity,
        number_of_bands, target_f0, frame_period, current_time, fs,
        world::kDefaultF0, &forward_real_fft, &inverse_real_fft,
        &minimum_phase, impulse_response);

    for (int i = current_position;
      i < MyMin(current_position + fft_size / 2, y_length - 1); ++i) {
      y[i] += impulse_response[i - current_position];
    }

    current_time += 1.0 / (f0[current_frame] ==
      0.0 ? world::kDefaultF0 : f0[current_frame]);
    current_frame = matlab_round(current_time / (frame_period / 1000.0));
    current_position = static_cast<int>(current_time * fs);
    if (current_frame >= f0_length) break;
  }

  DestroyMinimumPhaseAnalysis(&minimum_phase);
  DestroyInverseRealFFT(&inverse_real_fft);
  DestroyForwardRealFFT(&forward_real_fft);

  delete[] impulse_response;
}
