#include "ffluid.h"

static sim_data TmpLocal;
static sim_data RHSLocal[4];
const __float128 	one_sixth = 1.Q/6;
const __float128	one_half  = 0.5Q;
const __float128	two  =  2.0Q;
const __float128	one  =  1.0Q;

void ffluid_init_runge_kutta_4() {
  ffluid_data_init_copy(&DataCurr, &TmpLocal);
  ffluid_data_init_copy(&DataCurr, &RHSLocal[0]);
  ffluid_data_init_copy(&DataCurr, &RHSLocal[1]);
  ffluid_data_init_copy(&DataCurr, &RHSLocal[2]);
  ffluid_data_init_copy(&DataCurr, &RHSLocal[3]);
  ffluid_math_init_equations();
}

void ffluid_runge_kutta_4(data_ptr in) {
  EvolveConfig.cur_step++;
  in->time += EvolveConfig.dt;

  ffluid_data_copy (in, &TmpLocal);
  ffluid_call_rhs(&TmpLocal, &RHSLocal[0]);
  ffluid_data_fma(one_half*EvolveConfig.dt, &RHSLocal[0], in, &TmpLocal);
  ffluid_call_rhs(&TmpLocal, &RHSLocal[1]);
  ffluid_data_fma(one_half*EvolveConfig.dt, &RHSLocal[1], in, &TmpLocal);
  ffluid_call_rhs(&TmpLocal, &RHSLocal[2]);
  ffluid_data_fma(EvolveConfig.dt, &RHSLocal[2], in, &TmpLocal);
  ffluid_call_rhs(&TmpLocal, &RHSLocal[3]);

  // TmpLocal <- k3 + 2*k2 + 2*k1 + k0
  ffluid_data_fma(two, &RHSLocal[2], &RHSLocal[3], &TmpLocal);
  ffluid_data_fma(two, &RHSLocal[1], &TmpLocal, &TmpLocal);
  ffluid_data_fma(one, &RHSLocal[0], &TmpLocal, &TmpLocal);

  // out <- in + dt/6 * TmpLocal
  ffluid_data_fma(one_sixth*EvolveConfig.dt, &TmpLocal, in, in);
}

