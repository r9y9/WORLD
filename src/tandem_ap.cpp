//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise. All Rights Reserved.
// Author: morise [at] fc.ritsumei.ac.jp (Masanori Morise)
//
// Aperiodicity based on TANDEM-STRAIGHT.
// This function would be changed in near future.
//-----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include "./tandem_ap.h"
#include "./matlabfunctions.h"
#include "./constant_numbers.h"

const double kNormalCutoff = 600.0;

// Names of these variables are copied by the source code of Matlab.
// The developer does not know the meaning of these names.
typedef struct {
  int segmentLength;
  int nMargin;
  double **w;
  double *wsqrt;
  double **H;
  double **Hw;
  double **R;
  double **invR;

  double *Hwx;
  double *a;
  double *Ha;
  double *wxHa;
  double *wx;
} InternalParameters;

namespace {

//-----------------------------------------------------------------------------
// SetInternalParameters() allocates the memory to the struct.
//-----------------------------------------------------------------------------
void SetInternalParameters(int segment_length, int n_margin,
    InternalParameters *internal_parameters) {
  internal_parameters->segmentLength = segment_length;
  internal_parameters->nMargin = n_margin;
  internal_parameters->w = new double *[segment_length];
  for (int i = 0; i < segment_length; ++i)
    internal_parameters->w[i] = new double[segment_length];
  internal_parameters->wsqrt = new double[segment_length];
  internal_parameters->H = new double *[segment_length];
  for (int i = 0; i < segment_length; ++i)
    internal_parameters->H[i] = new double[n_margin * 2];
  internal_parameters->Hw = new double *[n_margin * 2];
  for (int i = 0; i < n_margin * 2; ++i)
    internal_parameters->Hw[i] = new double[segment_length];
  internal_parameters->R = new double *[n_margin * 2];
  for (int i = 0; i < n_margin * 2; ++i)
    internal_parameters->R[i] = new double[n_margin * 2];
  internal_parameters->invR = new double *[n_margin * 2];
  for (int i = 0; i < n_margin * 2; ++i)
    internal_parameters->invR[i] = new double[n_margin * 2];
  internal_parameters->Hwx = new double[n_margin * 2];
  internal_parameters->a = new double[n_margin * 2];
  internal_parameters->Ha = new double[segment_length];
  internal_parameters->wxHa = new double[segment_length];
  internal_parameters->wx = new double[segment_length];
}

//-----------------------------------------------------------------------------
// SetInternalParameters() frees the memory of the struct.
//-----------------------------------------------------------------------------
void DestroyInternalParameters(InternalParameters* internal_parameters) {
  delete[] internal_parameters->wsqrt;
  delete[] internal_parameters->wx;
  delete[] internal_parameters->wxHa;
  delete[] internal_parameters->Ha;
  delete[] internal_parameters->a;
  delete[] internal_parameters->Hwx;
  for (int i = 0; i < internal_parameters->nMargin * 2; ++i)
    delete[] internal_parameters->invR[i];
  delete[] internal_parameters->invR;
  for (int i = 0; i < internal_parameters->nMargin * 2; ++i)
    delete[] internal_parameters->R[i];
  delete[] internal_parameters->R;
  for (int i = 0; i < internal_parameters->nMargin * 2; ++i)
    delete[] internal_parameters->Hw[i];
  delete[] internal_parameters->Hw;
  for (int i = 0; i < internal_parameters->segmentLength; ++i)
    delete[] internal_parameters->H[i];
  delete[] internal_parameters->H;
  for (int i = 0; i < internal_parameters->segmentLength; ++i)
    delete[] internal_parameters->w[i];
  delete[] internal_parameters->w;
}

//-----------------------------------------------------------------------------
// Get*() calculates each parameter. The names do not followe the style guide.
// These names are refered by the article. To avoid the confusion,
// we employed the original names.
//-----------------------------------------------------------------------------
void GetH(double *x, int x_length, int segment_length, int index_bias,
    int current_position_in_sample, int t0_in_samples, double **H) {
  int index;
  for (int i = -1; i < 2; ++i) {
    for (int j = 0; j < segment_length; ++j) {
      index = std::max(0, std::min(x_length - 1,
        i + current_position_in_sample-index_bias - t0_in_samples + j));
      H[j][i + 1] = x[index];
      index = std::max(0, std::min(x_length - 1,
        i + current_position_in_sample - index_bias + t0_in_samples + j));
      H[j][i + 4] = x[index];
    }
  }
}

void GetHw(double **H, int segment_length, int n_margin2, double **w,
    double **Hw) {
  double tmp;
  for (int i = 0; i < n_margin2; ++i) {
    for (int j = 0; j < segment_length; ++j) {
      tmp = 0.0;
      for (int k = 0; k < segment_length; ++k) tmp += H[k][i] * w[k][j];
      Hw[i][j] = tmp;
    }
  }
}

void GetR(double **Hw, int n_margin2, int segment_length, double **H,
    double **R) {
  double tmp;
  for (int i = 0; i < n_margin2; ++i) {
    for (int j = 0; j < n_margin2; ++j) {
      tmp = 0.0;
      for (int k = 0; k < segment_length; ++k) tmp += Hw[i][k] * H[k][j];
      R[i][j] = tmp;
    }
  }
}

void GetHwx(double **Hw, int n_margin2, int segment_length, double *x,
    int origin, double *Hwx) {
  double tmp;
  for (int i = 0; i < n_margin2; ++i) {
    tmp = 0.0;
    for (int j = 0; j < segment_length; ++j) tmp += Hw[i][j]*x[origin+j];
    Hwx[i] = tmp;
  }
}

void Geta(double **invR, int n_margin2, double *Hwx, double *a) {
  double tmp;
  for (int i = 0; i < n_margin2; ++i) {
    tmp = 0.0;
    for (int j = 0; j < n_margin2; ++j) tmp += invR[i][j]*Hwx[j];
    a[i] = tmp;
  }
}

void GetHa(double **H, int segment_length, int n_margin2, double *a,
    double *Ha) {
  double tmp;
  for (int i = 0; i < segment_length; ++i) {
    tmp = 0.0;
    for (int j = 0; j < n_margin2; ++j) tmp += H[i][j]*a[j];
    Ha[i] = tmp;
  }
}

void GetW(int segment_length, double **w) {
  for (int i = 0; i < segment_length; ++i)
    for (int j = 0; j < segment_length; ++j) w[i][j] = 0.0;

  for (int i = 0; i < (segment_length - 1) / 2; ++i) {
    w[i][i] = 0.5 - 0.5 * cos((i + 1.0) /
      (segment_length + 1.0) * 2.0 * world::kPi);
    w[segment_length - i - 1][segment_length - i - 1] = w[i][i];
  }
  w[(segment_length - 1) / 2][(segment_length - 1) / 2] = 1.0;
}

double GetStdwxHa(double *wsqrt, int segment_length, double *x, int origin,
    double *Ha, double *wxHa) {
  for (int i = 0; i < segment_length; ++i)
    wxHa[i] = wsqrt[i] * (x[i + origin] - Ha[i]);
  return matlab_std(wxHa, segment_length);
}

double GetStdwx(double *wsqrt, int segment_length, double *x, int origin,
    double *wx) {
  for (int i = 0; i < segment_length; ++i) wx[i] = wsqrt[i] * x[i + origin];
  return matlab_std(wx, segment_length);
}

//-----------------------------------------------------------------------------
// f0PredictionResidualFixSegmentW() calculates the aperiodicity in a frequency
// band.
// This function is only used in BandwiseAperiodicity().
//-----------------------------------------------------------------------------
void f0PredictionResidualFixSegmentW(double *x, int x_length, double fs,
    double target_f0, double *temporalPositions, double *vuv, int f0_length,
    double initial_time, int duration_ms, int current_band,
    double **aperiodicity) {
  const int kNMargin = 3;
  int segment_length = matlab_round(fs * duration_ms / 2000.0) * 2 + 1;

  InternalParameters internal_parameters = {0};
  SetInternalParameters(segment_length, kNMargin, &internal_parameters);

  GetW(segment_length, internal_parameters.w);
  for (int i = 0; i < segment_length; ++i)
    internal_parameters.wsqrt[i] = sqrt(internal_parameters.w[i][i]);

  int t0_in_samples = matlab_round(fs / target_f0);
  int index_bias = matlab_round(fs / target_f0 / 2.0);

  int current_position_in_sample;
  int origin;
  for (int i = 0; i < f0_length; ++i) {
    current_position_in_sample =
      matlab_round(-initial_time + temporalPositions[i] * fs) + 1;
    if (vuv[i] != 0.0) {
      origin = std::max(0, std::min(x_length - 1,
        current_position_in_sample - index_bias));
      GetH(x, x_length, segment_length, index_bias,
          current_position_in_sample, t0_in_samples, internal_parameters.H);
      GetHw(internal_parameters.H, segment_length, kNMargin * 2,
          internal_parameters.w, internal_parameters.Hw);
      GetR(internal_parameters.Hw, kNMargin * 2, segment_length,
          internal_parameters.H, internal_parameters.R);
      GetHwx(internal_parameters.Hw, kNMargin * 2, segment_length,
          x, origin, internal_parameters.Hwx);
      inv(internal_parameters.R, kNMargin * 2, internal_parameters.invR);
      Geta(internal_parameters.invR, kNMargin * 2, internal_parameters.Hwx,
          internal_parameters.a);
      GetHa(internal_parameters.H, segment_length, kNMargin * 2,
          internal_parameters.a, internal_parameters.Ha);
      aperiodicity[i][current_band] = GetStdwxHa(internal_parameters.wsqrt,
          segment_length, x, origin, internal_parameters.Ha,
          internal_parameters.wxHa) / GetStdwx(internal_parameters.wsqrt,
          segment_length, x, origin, internal_parameters.wx);
    } else {  // Aperiodicity does not use if the speech is unvoiced.
      aperiodicity[i][current_band] = 0.0;
    }
  }
  DestroyInternalParameters(&internal_parameters);
}

//-----------------------------------------------------------------------------
// GetQMFpairOfFilters() sets the coefficients of QM filter (hHP:41, hLP:37)
// Although this function requires fs as the input, the result does not depend
// on it.
//-----------------------------------------------------------------------------
void GetQMFpairOfFilters(int fs, double *hHP, double *hLP) {
  // hHP
  hHP[0]  =  0.00041447996898231424;
  hHP[1]  =  0.00078125051417292477;
  hHP[2]  = -0.0010917236836275842;
  hHP[3]  = -0.0019867925675967589;
  hHP[4]  =  0.0020903896961562292;
  hHP[5]  =  0.0040940570272849346;
  hHP[6]  = -0.0034025808529816698;
  hHP[7]  = -0.0074961541272056016;
  hHP[8]  =  0.0049722633399330637;
  hHP[9]  =  0.012738791249119802;
  hHP[10] = -0.0066960326895749113;
  hHP[11] = -0.020694051570247052;
  hHP[12] =  0.0084324365650413451;
  hHP[13] =  0.033074383758700532;
  hHP[14] = -0.010018936738799522;
  hHP[15] = -0.054231361405808247;
  hHP[16] =  0.011293988915051487;
  hHP[17] =  0.10020081367388213;
  hHP[18] = -0.012120546202484579;
  hHP[19] = -0.31630021039095702;
  hHP[20] =  0.51240682580627639;
  hHP[21] = -0.31630021039095702;
  hHP[22] = -0.012120546202484579;
  hHP[23] =  0.10020081367388213;
  hHP[24] =  0.011293988915051487;
  hHP[25] = -0.054231361405808247;
  hHP[26] = -0.010018936738799522;
  hHP[27] =  0.033074383758700532;
  hHP[28] =  0.0084324365650413451;
  hHP[29] = -0.020694051570247052;
  hHP[30] = -0.0066960326895749113;
  hHP[31] =  0.012738791249119802;
  hHP[32] =  0.0049722633399330637;
  hHP[33] = -0.0074961541272056016;
  hHP[34] = -0.0034025808529816698;
  hHP[35] =  0.0040940570272849346;
  hHP[36] =  0.0020903896961562292;
  hHP[37] = -0.0019867925675967589;
  hHP[38] = -0.0010917236836275842;
  hHP[39] =  0.00078125051417292477;
  hHP[40] =  0.00041447996898231424;

  // hLP
  hLP[0]  = -0.00065488170077483048;
  hLP[1]  =  0.00007561994958159384;
  hLP[2]  =  0.0020408456937895227;
  hLP[3]  = -0.00074680535322030437;
  hLP[4]  = -0.0043502235688264931;
  hLP[5]  =  0.0025966428382642732;
  hLP[6]  =  0.0076396022827566962;
  hLP[7]  = -0.0064904118901497852;
  hLP[8]  = -0.011765804538954506;
  hLP[9]  =  0.013649908479276255;
  hLP[10] =  0.01636866479016021;
  hLP[11] = -0.026075976030529347;
  hLP[12] = -0.020910294856659444;
  hLP[13] =  0.048260725032316647;
  hLP[14] =  0.024767846611048111;
  hLP[15] = -0.096178467583360641;
  hLP[16] = -0.027359756709866623;
  hLP[17] =  0.31488052161630042;
  hLP[18] =  0.52827343594055032;
  hLP[19] =  0.31488052161630042;
  hLP[20] = -0.027359756709866623;
  hLP[21] = -0.096178467583360641;
  hLP[22] =  0.024767846611048111;
  hLP[23] =  0.048260725032316647;
  hLP[24] = -0.020910294856659444;
  hLP[25] = -0.026075976030529347;
  hLP[26] =  0.01636866479016021;
  hLP[27] =  0.013649908479276255;
  hLP[28] = -0.011765804538954506;
  hLP[29] = -0.0064904118901497852;
  hLP[30] =  0.0076396022827566962;
  hLP[31] =  0.0025966428382642732;
  hLP[32] = -0.0043502235688264931;
  hLP[33] = -0.00074680535322030437;
  hLP[34] =  0.0020408456937895227;
  hLP[35] =  0.00007561994958159384;
  hLP[36] = -0.00065488170077483048;
}

//-----------------------------------------------------------------------------
// GetSignalsForAperiodicity() calculates the signals used to calculate the
// aperiodicity. low_signal, high_signal and downsampled_high_signal are
// calculated in this function.
// This function is only used in BandwiseAperiodicity()
//-----------------------------------------------------------------------------
void GetSignalsForAperiodicity(int fft_size, double *whole_signal,
    int filtered_signal_length, double *hHP, double *hLP,
    double *low_signal, double *high_signal, double *downsampled_high_signal) {
  ForwardRealFFT forward_real_fft = {0};
  InverseRealFFT inverse_real_fft = {0};
  InitializeForwardRealFFT(fft_size, &forward_real_fft);
  InitializeInverseRealFFT(fft_size, &inverse_real_fft);
  fast_fftfilt(whole_signal, filtered_signal_length, hHP, 41,
      fft_size, &forward_real_fft, &inverse_real_fft, high_signal);
  fast_fftfilt(whole_signal, filtered_signal_length, hLP, 37,
      fft_size, &forward_real_fft, &inverse_real_fft, low_signal);
  DestroyForwardRealFFT(&forward_real_fft);
  DestroyInverseRealFFT(&inverse_real_fft);
  for (int j = 0; j < filtered_signal_length; j += 2)
    downsampled_high_signal[j / 2] = high_signal[j];
}

//-----------------------------------------------------------------------------
// UpdateWholeSignal() updates the whole_signal.
// This function is only used in BandwiseAperiodicity().
//-----------------------------------------------------------------------------
inline int UpdateWholeSignal(int filtered_signal_length, int fft_size,
    double *low_signal, double *whole_signal) {
  for (int i = 0; i < filtered_signal_length; i += 2)
    whole_signal[i / 2] = low_signal[i];
  for (int i = matlab_round(filtered_signal_length / 2.0); i < fft_size; i++)
    whole_signal[i] = 0.0;
  return matlab_round(filtered_signal_length / 2.0) + 82;
}

//-----------------------------------------------------------------------------
// BandwiseAperiodicity() calculates the aperiodicity in each frequency band.
//-----------------------------------------------------------------------------
void BandwiseAperiodicity(double *x, int x_length, int fs, double *f0,
    double *vuv, int f0_length, double *stretched_locations,
    int window_length_ms, double **aperiodicity) {
  double hHP[41], hLP[37];
  GetQMFpairOfFilters(fs, hHP, hLP);

  int number_of_bands =
    static_cast<int>(log(fs / kNormalCutoff) / world::kLog2);
  double *cutoff_list = new double[number_of_bands];

  for (int i = 0; i < number_of_bands; ++i)
    cutoff_list[i] = fs / pow(2.0, i + 2.0);

  // 82 = 41 (length of hHP) * 2
  int fft_size = GetSuitableFFTSize(x_length + 82);

  double *whole_signal = new double[fft_size];
  double *high_signal = new double[fft_size];
  double *low_signal = new double[fft_size];
  double *downsampled_high_signal = new double[fft_size];

  int filtered_signal_length = x_length + 82;

  for (int i = 0; i < x_length; ++i) whole_signal[i] = x[i];
  for (int i = x_length; i < fft_size; ++i) whole_signal[i] = 0.0;

  double tmp_fs;
  for (int i = 0; i < number_of_bands - 1; ++i) {
    tmp_fs = cutoff_list[i] * 2.0;
    GetSignalsForAperiodicity(fft_size, whole_signal, filtered_signal_length,
        hHP, hLP, low_signal, high_signal, downsampled_high_signal);

    f0PredictionResidualFixSegmentW(downsampled_high_signal,
        matlab_round(filtered_signal_length / 2.0), tmp_fs, f0[0],
        stretched_locations, vuv, f0_length, 41.0 / 2.0 / tmp_fs,
        window_length_ms, number_of_bands - i - 1, aperiodicity);

    // subband separation
    filtered_signal_length = UpdateWholeSignal(filtered_signal_length,
        fft_size, low_signal, whole_signal);
    // update the fft size
    fft_size = GetSuitableFFTSize(filtered_signal_length);
  }

  filtered_signal_length = (filtered_signal_length - 82) * 2;
  f0PredictionResidualFixSegmentW(whole_signal,
      matlab_round(filtered_signal_length / 2.0), tmp_fs, f0[0],
      stretched_locations, vuv, f0_length, 41.0 / 2.0 / tmp_fs,
      window_length_ms, 0, aperiodicity);

  delete[] downsampled_high_signal;
  delete[] low_signal;
  delete[] high_signal;
  delete[] whole_signal;
  delete[] cutoff_list;
}

//-----------------------------------------------------------------------------
// GetInterpolatedSignal() carries out the up sampling (target is 4 * fs)
//-----------------------------------------------------------------------------
void GetInterpolatedSignal(double *x, int x_length, double *interpolated_x) {
  interpolated_x[0] = x[0] * 0.14644660940672621;
  interpolated_x[1] = x[0] * 0.49999999999999994;
  interpolated_x[2] = x[0] * 0.85355339059327373;
  for (int i = 0; i < x_length - 1; ++i) {
    interpolated_x[i * 4 + 3] = x[i];
    interpolated_x[i * 4 + 4] = x[i] * 0.85355339059327373 +
      x[i + 1] * 0.14644660940672621;
    interpolated_x[i * 4 + 5] = x[i] * 0.49999999999999994 +
      x[i + 1] * 0.49999999999999994;
    interpolated_x[i * 4 + 6] = x[i] * 0.14644660940672621 +
      x[i + 1] * 0.85355339059327373;
  }
  interpolated_x[(x_length - 1) * 4 + 3] = x[x_length - 1];
  interpolated_x[(x_length - 1) * 4 + 4] =
    x[x_length - 1] * 0.85355339059327373;
  interpolated_x[(x_length - 1) * 4 + 5] =
    x[x_length - 1] * 0.49999999999999994;
  interpolated_x[(x_length - 1) * 4 + 6] =
    x[x_length - 1] * 0.14644660940672621;
  interpolated_x[(x_length - 1) * 4 + 7] =
    interpolated_x[(x_length - 1) * 4 + 8] =
    interpolated_x[(x_length - 1) * 4 + 9] = 0.0;
  return;
}

//-----------------------------------------------------------------------------
// GetNormalizedSignal() calculates the signal that the f0 contour is constant.
//-----------------------------------------------------------------------------
int GetNormalizedSignal(double *x, int x_length, int fs, double *f0,
    int f0_length, double frame_period, double target_f0,
    double **stretched_signal, double **stretched_locations) {
  int ix_length = x_length * 4 + 6;
//  int ix_length = x_length * 4;

  double *interpolated_x = new double[ix_length];
  GetInterpolatedSignal(x, x_length, interpolated_x);

  double *original_signal_time_axis = new double[ix_length];

  for (int i = 0; i < ix_length; ++i)
    original_signal_time_axis[i] = i / (fs * 4.0);

  double *base_f0 = new double[f0_length + 1];
  double *base_time_axis = new double[f0_length + 1];

  for (int i = 0; i < f0_length; ++i) {
    base_f0[i] = f0[i] == 0.0 ? target_f0 : f0[i];
    base_time_axis[i] = i * frame_period;
  }
  base_f0[f0_length] = target_f0;
  base_time_axis[f0_length] = f0_length * frame_period;

  double *interpolated_f0 = new double[ix_length];
  double *stretched_time_axis = new double[ix_length];
  interp1(base_time_axis, base_f0, f0_length + 1,
      original_signal_time_axis, ix_length, interpolated_f0);

  double tmp = target_f0 * fs * 4.0;
  stretched_time_axis[0] = interpolated_f0[0] / tmp;
  for (int i = 1; i < ix_length; ++i)
    stretched_time_axis[i] = stretched_time_axis[i - 1] +
    (interpolated_f0[i] / tmp);

  int stretched_signal_length =
    static_cast<int>(stretched_time_axis[ix_length - 1] * fs * 4.0) + 1;
  double *tmp_time_axis = new double[stretched_signal_length];
  double *stretched_signal4 = new double[stretched_signal_length];

  for (int i = 0; i < stretched_signal_length; ++i)
    tmp_time_axis[i] = i / (fs * 4.0);
  interp1(stretched_time_axis, interpolated_x, ix_length,
    tmp_time_axis, stretched_signal_length, stretched_signal4);

  *stretched_locations = new double[f0_length];
  interp1(original_signal_time_axis, stretched_time_axis, ix_length,
    base_time_axis, f0_length, *stretched_locations);

  // 17 is a safe guard.
  *stretched_signal = new double[stretched_signal_length / 4 + 17];
  decimate(stretched_signal4, stretched_signal_length, 4, *stretched_signal);

  delete[] stretched_signal4;
  delete[] tmp_time_axis;
  delete[] stretched_time_axis;
  delete[] base_f0;
  delete[] base_time_axis;
  delete[] interpolated_f0;
  delete[] original_signal_time_axis;
  delete[] interpolated_x;

  return 1 + stretched_signal_length / 4;
}

}  // namespace

int GetNumberOfBands(int fs) {
  return static_cast<int>(log(fs / kNormalCutoff) / world::kLog2);
}

double AperiodicityRatio(double *x, int x_length, int fs, double *f0,
    int f0_length, double frame_period, double **aperiodicity) {
  double max_f0 = 0.0;
  for (int i = 0; i < f0_length; ++i)
    max_f0 = max_f0 > f0[i] ? max_f0 : f0[i];

  const double kMinimumF0ForNormalization = 32.0;
  const double kMaximumF0ForNormalization = 200.0;
  double target_f0 = std::max(kMinimumF0ForNormalization,
      std::min(kMaximumF0ForNormalization, max_f0));

  // The number of two arraies are unknown.
  double *stretched_signal = NULL;
  double *stretched_locations = NULL;

  int normalized_signal_length = GetNormalizedSignal(x, x_length, fs, f0,
      f0_length, frame_period / 1000.0, target_f0, &stretched_signal,
      &stretched_locations);

  double *stretched_f0 = new double[f0_length];
  for (int i = 0; i < f0_length; ++i) stretched_f0[i] = target_f0;

  BandwiseAperiodicity(stretched_signal, normalized_signal_length, fs,
      stretched_f0, f0, f0_length, stretched_locations,
      matlab_round(2000.0 / target_f0), aperiodicity);

  delete[] stretched_f0;
  delete[] stretched_signal;
  delete[] stretched_locations;

  return target_f0;
}
