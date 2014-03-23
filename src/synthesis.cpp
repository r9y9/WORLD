//-----------------------------------------------------------------------------
// Copyright 2012-2013 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//
// Voice synthesis based on f0, spectrogram and spectrogram of
// excitation signal.
//-----------------------------------------------------------------------------
#include "./synthesis.h"

#include <math.h>
#include <stdlib.h>

#include "./common.h"
#include "./constantnumbers.h"
#include "./matlabfunctions.h"

namespace {

//-----------------------------------------------------------------------------
// GetOneFrameSegment() calculates a glottal vibration based on the spectral
// envelope and excitation signal.
// Caution:
//   minimum_phase and inverse_real_fft are allocated in advance. This is for
//   the rapid processing because set of FFT requires much computational cost.
//-----------------------------------------------------------------------------
void GetOneFrameSegment(double *f0, double **spectrogram,
    double **residual_spectrogram, int fft_size, int current_frame,
    MinimumPhaseAnalysis *minimum_phase, InverseRealFFT *inverse_real_fft,
    double *y) {
  for (int i = 0; i <= minimum_phase->fft_size / 2; ++i)
    minimum_phase->log_spectrum[i] =
    log(spectrogram[current_frame][i]) / 2.0;
  GetMinimumPhaseSpectrum(minimum_phase);

  inverse_real_fft->spectrum[0][0] =
    minimum_phase->minimum_phase_spectrum[0][0] *
    residual_spectrogram[current_frame][0];
  inverse_real_fft->spectrum[0][1] = 0;

  for (int i = 1; i < fft_size / 2; ++i) {
    inverse_real_fft->spectrum[i][0] =
      minimum_phase->minimum_phase_spectrum[i][0] *
      residual_spectrogram[current_frame][(i - 1) * 2 + 1] -
      minimum_phase->minimum_phase_spectrum[i][1] *
      residual_spectrogram[current_frame][i * 2];
    inverse_real_fft->spectrum[i][1] =
      minimum_phase->minimum_phase_spectrum[i][0] *
      residual_spectrogram[current_frame][i * 2] +
      minimum_phase->minimum_phase_spectrum[i][1] *
      residual_spectrogram[current_frame][(i - 1) * 2 + 1];
  }
  inverse_real_fft->spectrum[fft_size / 2][0] =
    minimum_phase->minimum_phase_spectrum[fft_size / 2][0] *
    residual_spectrogram[current_frame][fft_size - 1];
  inverse_real_fft->spectrum[fft_size / 2][1] = 0;
  fft_execute(inverse_real_fft->inverse_fft);

  for (int i = 0; i < fft_size; ++i)
    y[i] = inverse_real_fft->waveform[i] / fft_size;
}

}  // namespace

void Synthesis(double *f0, int f0_length, double **spectrogram,
    double **residual_spectrogram, int fft_size, double frame_period,
    int fs, int y_length, double *y) {
  double *impulse_response = new double[fft_size];

  MinimumPhaseAnalysis minimum_phase = {0};
  InitializeMinimumPhaseAnalysis(fft_size, &minimum_phase);
  InverseRealFFT inverse_real_fft = {0};
  InitializeInverseRealFFT(fft_size, &inverse_real_fft);

  double current_time = 0.0;
  int current_position = 0;
  int current_frame = 0;
  // Length used for the synthesis is unclear.
  const int kFrameLength = 3 * fft_size / 4;

  while (1) {
    GetOneFrameSegment(f0, spectrogram, residual_spectrogram, fft_size,
        current_frame, &minimum_phase, &inverse_real_fft, impulse_response);

    for (int i = current_position;
        i < MyMin(current_position + kFrameLength, y_length - 1); ++i)
      y[i] += impulse_response[i - current_position];

    // update
    current_time +=
      1.0 / (f0[current_frame] == 0.0 ? world::kDefaultF0 : f0[current_frame]);
    current_frame = matlab_round(current_time / (frame_period / 1000.0));
    current_position = static_cast<int>(current_time * fs);
    if (current_frame >= f0_length) break;
  }

  DestroyMinimumPhaseAnalysis(&minimum_phase);
  DestroyInverseRealFFT(&inverse_real_fft);
  delete[] impulse_response;
}
