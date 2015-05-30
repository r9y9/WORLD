//-----------------------------------------------------------------------------
// Copyright 2012-2015 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//
//-----------------------------------------------------------------------------
#ifndef WORLD_DIO_H_
#define WORLD_DIO_H_

#include "./dllexport.h"

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
// Struct for DIO
//-----------------------------------------------------------------------------
typedef struct {
  double f0_floor;
  double f0_ceil;
  double channels_in_octave;
  double frame_period;  // msec
  int speed;  // (1, 2, ..., 12)
} DioOption;

//-----------------------------------------------------------------------------
// DIO
// Input:
//   x          : Input signal
//   x_length   : Length of x
//   fs         : Sampling frequency
//   option     : Struct to order the parameter for DIO
// Output:
//   time_axis  : Temporal positions.
//   f0         : F0 contour.
//-----------------------------------------------------------------------------
DLLEXPORT void Dio(double *x, int x_length, int fs, const DioOption option,
  double *time_axis, double *f0);


// This function is basically same as Dio but different in argument type:
// `const DioOption` -> `const DioOption*`
// This function can be called by languages that doesn't support struct passing
// by value (e.g. Julia)
DLLEXPORT void DioByOptPtr(double *x, int x_length, int fs, const DioOption *option,
  double *time_axis, double *f0);

//-----------------------------------------------------------------------------
// InitializeDioOption allocates the memory to the struct and sets the
// default parameters.
// Output:
//   option   : Struct for the optional parameter.
//-----------------------------------------------------------------------------
DLLEXPORT void InitializeDioOption(DioOption *option);

//-----------------------------------------------------------------------------
// GetSamplesForDIO() calculates the number of samples required for Dio().
// Input:
//   fs             : Sampling frequency [Hz]
//   x_length       : Length of the input signal [Sample].
//   frame_period   : Frame shift [msec]
// Output:
//   The number of samples required to store the results of Dio()
//-----------------------------------------------------------------------------
DLLEXPORT int GetSamplesForDIO(int fs, int x_length, double frame_period);

#ifdef __cplusplus
}
#endif

#endif  // WORLD_DIO_H_
