//-----------------------------------------------------------------------------
// Copyright 2012-2014 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//-----------------------------------------------------------------------------
#ifndef WORLD_CHEAPTRICK_H_
#define WORLD_CHEAPTRICK_H_

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
// CheapTrick() calculates the spectrogram that consists of spectral envelopes
// estimated by CheapTrick.
// Input:
//   x            : Input signal
//   x_length     : Length of x
//   fs           : Sampling frequency
//   time_axis    : Time axis
//   f0           : F0 contour
// Output:
//   spectrogram  : Spectrogram estimated by STAR.
//-----------------------------------------------------------------------------
void CheapTrick(double *x, int x_length, int fs, double *time_axis, double *f0,
  int f0_length, double **spectrogram);

//-----------------------------------------------------------------------------
// GetFFTSizeForCheapTrick() calculates the FFT size based on the sampling
// frequency and the lower limit of f0 (It is defined in world.h).
// Input:
//   fs      : Sampling frequency
// Output:
//   FFT size
//-----------------------------------------------------------------------------
int GetFFTSizeForCheapTrick(int fs);

#ifdef __cplusplus
}
#endif

#endif  // WORLD_CHEAPTRICK_H_
