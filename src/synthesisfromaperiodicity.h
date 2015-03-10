//-----------------------------------------------------------------------------
// Copyright 2012-2014 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//-----------------------------------------------------------------------------
#ifndef WORLD_SYNTHESISFROMAPERIODICITY_H_
#define WORLD_SYNTHESISFROMAPERIODICITY_H_

#include "./dllexport.h"

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
// synthesis_ap() synthesize the voice based on f0, spectrogram and
// aperiodicity (not excitation signal).
// Input:
//   f0                   : f0 contour
//   f0_length            : Length of f0
//   spectrogram          : Spectrogram estimated by STAR
//   fft_size             : FFT size
//   aperiodicity         : Aperiodicity spectrogram based on TANDEM_AP
//   frame_period         : Temporal period used for the analysis
//   fs                   : Sampling frequency
//   y_length             : Length of the output signal (Memory of y has been
//                          allocated in advance)
// Output:
//   y                    : Calculated voice
//-----------------------------------------------------------------------------
DLLEXPORT void SynthesisFromAperiodicity(double *f0, int f0_length,
    double **spectrogram, double **aperiodicity, int fft_size,
    double frame_period, int fs, int y_length, double *y);

#ifdef __cplusplus
}
#endif

#endif  // WORLD_SYNTHESISFROMAPERIODICITY_H_
