//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise. All Rights Reserved.
// Author: morise [at] fc.ritsumei.ac.jp (Masanori Morise)
//
// common.cpp includes functions used in at least two files.
// (1) Common functions
// (2) FFT, IFFT and minimum phase analysis.
//
// In FFT analysis and minimum phase analysis,
// Functions "Initialize*()" allocate the mamory.
// Functions "Destroy*()" free the accolated memory.
// FFT size is used for initialization, and structs are used to keep the memory.
// Functions "GetMinimumPhaseSpectrum()" calculate minimum phase spectrum.
// Forward and inverse FFT do not have the function "Get*()",
// because forward FFT and inverse FFT can run in one step.
//
//-----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "./common.h"
#include "./constant_numbers.h"

//-----------------------------------------------------------------------------
// Fundamental functions

int GetSuitableFFTSize(int sample) {
  return static_cast<int>(pow(2.0,
    static_cast<int>(log(static_cast<double>(sample)) / world::kLog2) + 1.0));
}

//-----------------------------------------------------------------------------
// FFT, IFFT and minimum phase analysis
void InitializeForwardRealFFT(int fft_size, ForwardRealFFT *forward_real_fft) {
  forward_real_fft->fft_size = fft_size;
  forward_real_fft->waveform = new double[fft_size];
  forward_real_fft->spectrum = new fft_complex[fft_size];
  forward_real_fft->forward_fft = fft_plan_dft_r2c_1d(fft_size,
      forward_real_fft->waveform, forward_real_fft->spectrum, FFT_ESTIMATE);
}

void DestroyForwardRealFFT(ForwardRealFFT *forward_real_fft) {
  fft_destroy_plan(forward_real_fft->forward_fft);
  delete[] forward_real_fft->spectrum;
  delete[] forward_real_fft->waveform;
}

void InitializeInverseRealFFT(int fft_size, InverseRealFFT *inverse_real_fft) {
  inverse_real_fft->fft_size = fft_size;
  inverse_real_fft->waveform = new double[fft_size];
  inverse_real_fft->spectrum = new fft_complex[fft_size];
  inverse_real_fft->inverse_fft = fft_plan_dft_c2r_1d(fft_size,
      inverse_real_fft->spectrum, inverse_real_fft->waveform, FFT_ESTIMATE);
}

void DestroyInverseRealFFT(InverseRealFFT *inverse_real_fft) {
  fft_destroy_plan(inverse_real_fft->inverse_fft);
  delete[] inverse_real_fft->spectrum;
  delete[] inverse_real_fft->waveform;
}

void InitializeMinimumPhaseAnalysis(int fft_size,
                                    MinimumPhaseAnalysis *minimum_phase) {
  minimum_phase->fft_size = fft_size;
  minimum_phase->log_spectrum = new double[fft_size];
  minimum_phase->minimum_phase_spectrum = new fft_complex[fft_size];
  minimum_phase->cepstrum = new fft_complex[fft_size];
  minimum_phase->inverse_fft = fft_plan_dft_r2c_1d(fft_size,
      minimum_phase->log_spectrum, minimum_phase->cepstrum, FFT_ESTIMATE);
  minimum_phase->forward_fft = fft_plan_dft_1d(fft_size,
      minimum_phase->cepstrum, minimum_phase->minimum_phase_spectrum,
      FFT_FORWARD, FFT_ESTIMATE);
}

void GetMinimumPhaseSpectrum(MinimumPhaseAnalysis *minimum_phase) {
  // Mirroring
  for (int i = minimum_phase->fft_size / 2 + 1;
       i < minimum_phase->fft_size; ++i)
    minimum_phase->log_spectrum[i] =
    minimum_phase->log_spectrum[minimum_phase->fft_size - i];

  // This fft_plan carries out "forward" FFT.
  // To carriy out the Inverse FFT, the sign of imaginary part
  // is inverted after FFT.
  fft_execute(minimum_phase->inverse_fft);
  minimum_phase->cepstrum[0][1] *= -1.0;
  for (int i = 1; i < minimum_phase->fft_size / 2; ++i) {
    minimum_phase->cepstrum[i][0] *= 2.0;
    minimum_phase->cepstrum[i][1] *= -2.0;
  }
  minimum_phase->cepstrum[minimum_phase->fft_size / 2][1] *= -1.0;
  for (int i = minimum_phase->fft_size / 2 + 1;
       i < minimum_phase->fft_size; ++i) {
    minimum_phase->cepstrum[i][0] = 0.0;
    minimum_phase->cepstrum[i][1] = 0.0;
  }
  fft_execute(minimum_phase->forward_fft);

  // This FFT library does not keep the aliasing.
  // Since x is complex number, calculation of exp(x) is as following.
  double tmp;
  for (int i = 0; i <= minimum_phase->fft_size / 2; ++i) {
    tmp = exp(minimum_phase->minimum_phase_spectrum[i][0] /
      minimum_phase->fft_size);
    minimum_phase->minimum_phase_spectrum[i][0] = tmp *
      cos(minimum_phase->minimum_phase_spectrum[i][1] /
      minimum_phase->fft_size);
    minimum_phase->minimum_phase_spectrum[i][1] = tmp *
      sin(minimum_phase->minimum_phase_spectrum[i][1] /
      minimum_phase->fft_size);
  }
}

void DestroyMinimumPhaseAnalysis(MinimumPhaseAnalysis *minimum_phase) {
  fft_destroy_plan(minimum_phase->forward_fft);
  fft_destroy_plan(minimum_phase->inverse_fft);
  delete[] minimum_phase->cepstrum;
  delete[] minimum_phase->log_spectrum;
  delete[] minimum_phase->minimum_phase_spectrum;
}

