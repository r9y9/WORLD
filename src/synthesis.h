//-----------------------------------------------------------------------------
// Copyright 2012-2014 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//-----------------------------------------------------------------------------
#ifndef WORLD_SYNTHESIS_H_
#define WORLD_SYNTHESIS_H_

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
// Synthesis() synthesizes the voice from three parameters.
// Input:
//   f0                     : f0 contour
//   f0_length              : Length of f0
//   spectrogram            : Spectrogram (WORLD assumes spectrogram by Star())
//   residual_spectrogram   : Extracted spectrum of the excitation signal
//   fft_size               : FFT size used for Star() and Platinum()
//   frame_period           : Temporal inverval for Dio()
//   fs                     : Sampling frequency
//   y_length               : Length of Output (allocated in advance)
// Output:
//   y                      : Synthesized voice
//-----------------------------------------------------------------------------
void Synthesis(double *f0, int f0_length, double **spectrogram,
  double **residual_spectrogram, int fft_size, double frame_period, int fs,
  int y_length, double *y);

#ifdef __cplusplus
}
#endif

#endif  // WORLD_SYNTHESIS_H_
