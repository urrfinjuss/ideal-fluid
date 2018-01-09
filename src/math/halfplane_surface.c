#include "ffluid.h"

static sim_data SimLocal;
static aux_data	AuxLocal;
static fft_list FFTLocal; 

void ffluid_math_init_surface_halfplane() {
  ffluid_data_init_copy(&DataCurr, &SimLocal);
  ffluid_alloc_aux_array(&AuxLocal, 2, DataCurr.N);
  ffluid_alloc_fft_plans(&AuxLocal, &FFTLocal);
}

void ffluid_math_clear_surface_halfplane() {
  ffluid_dealloc_aux_array(&AuxLocal);
}


void ffluid_math_get_surface_variables_halfplane(data_ptr in, data_ptr out) {
  /* get surface potential \Phi(q) and z(q), but the zero modes are not set. */
  /* auxiliary variable S = z_u u_q is introduced */
  unsigned long 	N = in->N;
  long_double_t		S0 = 0.0L, T0 = 0.0L;
  long_complex_t	z0 = 0.L;
  
  /* direct inverse O(N log N) */
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.ph[0][j] = (cpowl(in->Q[j], -2) - 1.0L)*in->du[j]/N;
    AuxLocal.ph[1][j] = -1.0IL*in->V[j]*(AuxLocal.ph[0][j] + 1.0L)*in->du[j]/N;
  }
  /* can be improved with advanced FFTW interface */
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  memset(AuxLocal.ft[0]+N/2, 0, N/2);
  memset(AuxLocal.ft[1]+N/2, 0, N/2);
  AuxLocal.ft[0][0] = 0.0L;
  AuxLocal.ft[1][0] = 0.0L;
  for (unsigned long j = N/2-1; j > 0; j--) {
    /* AuxLocal.ft[0] stores Fourier modes of z_k, k < 0 */
    /* AuxLocal.ft[1] stores Fourier modes of Phi_k, k < 0 */ 
    AuxLocal.ft[0][j] = 1.0IL*AuxLocal.ft[0][j]/j;
    AuxLocal.ft[1][j] = 1.0IL*AuxLocal.ft[1][j]/j;
    S0 += j*creall(AuxLocal.ft[0][j]*conjl(AuxLocal.ft[0][j]));
    T0 += -creall(AuxLocal.ft[0][j]);
  }
  /* setting of the zero modes of Phi,Z goes here */
  memcpy(SimLocal.Q, AuxLocal.ft[0], N*sizeof(long_complex_t));
  memcpy(SimLocal.V, AuxLocal.ft[1], N*sizeof(long_complex_t));
  ffluid_math_set_zero_mode(&SimLocal, -0.5IL*S0, &z0);
  AuxLocal.ft[0][0] = T0 - creall(z0) + z0;
  /* can be improved with advanced FFTW interface */ 
  fftwl_execute(FFTLocal.bp[0]);
  fftwl_execute(FFTLocal.bp[1]);
  for (unsigned long j = 0; j < N; j++) {
    out->Q[j] = in->u[j] + AuxLocal.ph[0][j];
  }
  memcpy(out->V, AuxLocal.ph[1], N*sizeof(long_complex_t));
}


void ffluid_math_set_zero_mode(data_ptr in, long_complex_t S0, long_complex_t *out) {
  unsigned long		N = in->N;
  long_complex_t	w = cexpl(1.IL*in->q0);
  long_complex_t	tmp;
  long_double_t 	b = 0.5L*(1.L + powl(in->l, 2))/in->l;
  long_double_t 	xi = (1.L - powl(in->l, 2))/(1.L + powl(in->l, 2));

  //tmpc[2][0] = -0.5L*xi*conjl(w);
  AuxLocal.ph[0][0] = -0.5L*xi*conjl(w);
  /*
  for (long int j = 1; j < state.number_modes/2-1; j++) {
    tmpc[2][j] = tmpc[2][0]/(1.L - conjl(tmpc[2][0])*tmpc[2][j-1]); 
  }
  */
  for (unsigned long j = 1; j < N/2 - 1; j++) {
    AuxLocal.ph[0][j] = AuxLocal.ph[0][0]/(1.L - conjl(AuxLocal.ph[0][0])*AuxLocal.ph[0][j-1]);
  }
  //tmpc[3][0] = in[1]/b - 2.L*S0*conjl(tmpc[2][0]);
  AuxLocal.ph[1][0] = in->Q[1]/b - 2.0L*S0*conjl(AuxLocal.ph[0][0]);
  /*
  for (long int j = 1; j < state.number_modes/2-1; j++) {
    tmpc[3][j] = in[j+1]/b - conjl(tmpc[2][0])*tmpc[3][j-1];
    tmpc[3][j] = tmpc[3][j]/(1.L - conjl(tmpc[2][0])*tmpc[2][j-1]);
  } 
  */
  for (unsigned long j = 1; j < N/2 - 1; j++) {
    AuxLocal.ph[1][j] = in->Q[j+1]/b - conjl(AuxLocal.ph[0][0])*AuxLocal.ph[1][j-1];
    AuxLocal.ph[1][j] = AuxLocal.ph[1][j]/(1.0L - conjl(AuxLocal.ph[0][0])*AuxLocal.ph[0][j-1]);
  }

  //tmpc[4][state.number_modes/2-2] = tmpc[3][state.number_modes/2-2];
  tmp = AuxLocal.ph[1][N/2-2];
  /*
  for (long int j = state.number_modes/2-3; j > -1; j--) {
    tmpc[4][j] = tmpc[3][j]-tmpc[2][j]*tmpc[4][j+1];
  }
  */
  for (long int j = N/2-3; j > -1; j--) {
    tmp = AuxLocal.ph[1][j]-AuxLocal.ph[0][j]*tmp;
  }
  *out = b*(tmp*AuxLocal.ph[0][0] + 1.0L*S0);
  //*out = b*(tmpc[4][0]*tmpc[2][0] + 1.L*S0);
}



