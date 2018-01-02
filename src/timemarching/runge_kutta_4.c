#include "ffluid.h"

data_ptr tmp;
data_ptr *rhs;
const __float128 	one_sixth = 1.Q/6;
const long_double_t	one_half  = (long_double_t) 0.5Q;
const long_double_t	two  = (long_double_t) 2.0Q;
const long_double_t	one  = (long_double_t) 1.0Q;

void ffluid_runge_kutta_4(data_ptr in) {
  EvolveConfig.cur_step++;
  in->time += EvolveConfig.dt;

  ffluid_data_copy (in, tmp);
  ffluid_call_rhs(tmp, rhs[0]);
  ffluid_data_fma(one_half*EvolveConfig.dt, rhs[0], in, tmp);
  ffluid_call_rhs(tmp, rhs[1]);
  ffluid_data_fma(one_half*EvolveConfig.dt, rhs[1], in, tmp);
  ffluid_call_rhs(tmp, rhs[2]);
  ffluid_data_fma(EvolveConfig.dt, rhs[2], in, tmp);
  ffluid_call_rhs(tmp, rhs[3]);

  // tmp <- k3 + 2*k2 + 2*k1 + k0
  ffluid_data_fma(two, rhs[2], rhs[3], tmp);
  ffluid_data_fma(two, rhs[1], tmp, tmp);
  ffluid_data_fma(one, rhs[0], tmp, tmp);

  // out <- in + dt/6 * tmp
  ffluid_data_fma(one_sixth*EvolveConfig.dt, tmp, in, in);
}

