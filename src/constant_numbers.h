//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise. All Rights Reserved.
// Author: morise [at] fc.ritsumei.ac.jp (Masanori Morise)
//
// This header file only defines constant numbers used for several function.
//-----------------------------------------------------------------------------

#ifndef WORLD_CONSTANT_NUMBERS_H_
#define WORLD_CONSTANT_NUMBERS_H_

#ifdef __cplusplus
extern "C" {
#endif

namespace world {
  const double kPi = 3.1415926535897932384;
  const double kMySafeGuardMinimum = 0.000000000001;
  const double kFloorF0 = 71.0;
  const double kDefaultF0 = 150.0;
  const double kLog2 = 0.69314718055994529;
  // Maximum standard deviation not to be selected as a best f0.
  const double kMaximumValue = 100000.0;
// Note to me (fs: 44100)
// 71 Hz is the limit to maintain the FFT size at 2048.
// If we use 70 Hz as FLOOR_F0, the FFT size of 4096 is required.

}  // namespace world

#ifdef __cplusplus
}
#endif

#endif
