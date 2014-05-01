//-----------------------------------------------------------------------------
// Copyright 2012-2013 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//-----------------------------------------------------------------------------
#ifndef WORLD_TANDEM_AP_H_
#define WORLD_TANDEM_AP_H_

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
// GetNumberOfBands() calculate the number of bands for aperiodicity.
// Input:
//   fs       : Sampling frequency
// Output:
//   The number of bands required for the calculation
//-----------------------------------------------------------------------------
int GetNumberOfBands(int fs);

//-----------------------------------------------------------------------------
// The latest version of aperiodicity estimation in TANDEM-STRAIGHT.
// This function skipped several complex processes.
// Input:
//   x            : Input signal
//   x_length     : Length of x
//   f0           : f0 contour
//   f0_length    : Length of f0
//   frame_period : Time interval for analysis
// Output:
//   aperiodicity : Estimated aperiodicity
//   Value used for the aperiodicity estimation. This value is used for
//   the synthesis.
//-----------------------------------------------------------------------------
double AperiodicityRatio(double *x, int x_length, int fs, double *f0,
  int f0_length, double frame_period, double **aperiodicity);

#ifdef __cplusplus
}
#endif

#endif  // WORLD_TANDEM_AP_H_
