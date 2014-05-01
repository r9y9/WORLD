//-----------------------------------------------------------------------------
// Copyright 2012-2014 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//-----------------------------------------------------------------------------
#ifndef WORLD_STONEMASK_H_
#define WORLD_STONEMASK_H_

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
// StoneMask() refines the estimated F0 by Dio()
// Input:
//   x                      : Input signal
//   x_length               : Length of the input signal
//   fs                     : Sampling frequency
//   time_axis              : Temporal information
//   f0                     : f0 contour
//   f0_length              : Length of f0
// Output:
//   refined_f0             : Refined F0
//-----------------------------------------------------------------------------
void StoneMask(double *x, int x_length, int fs, double *time_axis, double *f0,
    int f0_length, double *refined_f0);

#ifdef __cplusplus
}
#endif

#endif  // WORLD_STONEMASK_H_
