//-----------------------------------------------------------------------------
// Copyright 2012-2013 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//-----------------------------------------------------------------------------
#ifndef WORLD_SYNTHESIS_AP_H_
#define WORLD_SYNTHESIS_AP_H_

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
//   number_of_bands      : Number of frequency bands used for TANDEM_AP
//   target_f0            : Only a parameter in TANDEM_AP
//   frame_period         : Temporal period used for the analysis
//   fs                   : Sampling frequency
//   y_length             : Length of the output signal (Memory of y has been
//                          allocated in advance)
// Output:
//   y                    : Calculated glottal pulse
//-----------------------------------------------------------------------------
void SynthesisFromAperiodicity(double *f0, int f0_length, double **spectrogram,
  int fft_size, double **aperiodicity, int number_of_bands, double target_f0,
  double frame_period, int fs, int y_length, double *synthesisOut);

#ifdef __cplusplus
}
#endif

#endif  // WORLD_SYNTHESIS_AP_H_
