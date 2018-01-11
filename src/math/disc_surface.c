#include "ffluid.h"

static sim_data SimLocal;
static aux_data AuxLocal;
static fft_list FFTLocal;

void ffluid_math_init_surface() {
  ffluid_data_init_copy(&DataCurr, &SimLocal);
  ffluid_alloc_aux_array(&AuxLocal, 2, DataCurr.N);
  ffluid_alloc_fft_plans(&AuxLocal, &FFTLocal);
}

void ffluid_math_get_surface_variables(data_ptr in, data_ptr out) { 
}
