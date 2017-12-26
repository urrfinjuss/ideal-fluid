#include "ffluid.h"

data_ptr tmp;
data_ptr *rhs;
const double one_sixth = 1./6;

void runge_kutta_4(data_ptr in, const double step_size) {
  ffluid_data_copy (in, tmp);
  ffluid_call_rhs(tmp, rhs[0]);
  ffluid_data_fma(0.5*step_size, rhs[0], in, tmp);
  ffluid_call_rhs(tmp, rhs[1]);
  ffluid_data_fma(0.5*step_size, rhs[1], in, tmp);
  ffluid_call_rhs(tmp, rhs[2]);
  ffluid_data_fma(1.0*step_size, rhs[2], in, tmp);
  ffluid_call_rhs(tmp, rhs[3]);

  // tmp <- k3 + 2*k2 + 2*k1 + k0
  ffluid_data_fma(2.0, rhs[2], rhs[3], tmp);
  ffluid_data_fma(2.0, rhs[1], tmp, tmp);
  ffluid_data_fma(1.0, rhs[0], tmp, tmp);

  // out <- in + dt/6 * tmp
  ffluid_data_fma(one_sixth*step_size, tmp, in, in);
}

