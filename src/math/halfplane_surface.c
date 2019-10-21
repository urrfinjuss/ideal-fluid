#include "ffluid.h"

static sim_data SimLocal;
static aux_data	AuxLocal;
static fft_list FFTLocal; 

void ffluid_math_init_surface() {
  ffluid_data_init_copy(&DataCurr, &SimLocal);
  ffluid_alloc_aux_array(&AuxLocal, 2, DataCurr.N);
  ffluid_alloc_fft_plans(&AuxLocal, &FFTLocal);
}

void ffluid_math_clear_surface() {
  ffluid_dealloc_aux_array(&AuxLocal);
}

void ffluid_math_get_r0(data_ptr in) {
}

void ffluid_math_get_surface_variables_RV(data_ptr in, data_ptr out) {
}
void ffluid_math_get_surface_variables(data_ptr in, data_ptr out) {
  /* auxiliary variable S = z_u u_q is introduced */
  unsigned long 	N = in->N;
  long_double_t		S0 = 0.0L, T0 = 0.0L;
  long_complex_t	z0 = 0.L;
  
  /* direct inverse O(N log N) */
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[0][j] = (cpowl(in->R[j], -1) - 1.0L)*in->du[j]/N;
    AuxLocal.X[1][j] = -1.0IL*in->V[j]*(AuxLocal.X[0][j] + 1.0L)*in->du[j]/N;
  }
  /* can be improved with advanced FFTW interface */
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  memset(AuxLocal.Y[0]+N/2, 0, N/2);
  memset(AuxLocal.Y[1]+N/2, 0, N/2);
  AuxLocal.Y[0][0] = 0.0L;
  AuxLocal.Y[1][0] = 0.0L;
  for (unsigned long j = N/2-1; j > 0; j--) {
    /* AuxLocal.ft[0] stores Fourier modes of z_k, k < 0 */
    /* AuxLocal.ft[1] stores Fourier modes of Phi_k, k < 0 */ 
    AuxLocal.Y[0][j] = 1.0IL*AuxLocal.Y[0][j]/j;
    AuxLocal.Y[1][j] = 1.0IL*AuxLocal.Y[1][j]/j;
    S0 += j*creall(AuxLocal.Y[0][j]*conjl(AuxLocal.Y[0][j]));
    T0 += -creall(AuxLocal.Y[0][j]);
  }
  /* setting of the zero modes of Phi,Z goes here */
  memcpy(SimLocal.R, AuxLocal.Y[0], N*sizeof(long_complex_t));
  memcpy(SimLocal.V, AuxLocal.Y[1], N*sizeof(long_complex_t));
  ffluid_math_set_zero_mode(&SimLocal, -0.5IL*S0, &z0);
  AuxLocal.Y[0][0] = T0 - creall(z0) + z0;
  /* can be improved with advanced FFTW interface */ 
  fftwl_execute(FFTLocal.bp[0]);
  fftwl_execute(FFTLocal.bp[1]);
  for (unsigned long j = 0; j < N; j++) {
    out->R[j] = in->u[j] + AuxLocal.X[0][j];
  }
  memcpy(out->V, AuxLocal.X[1], N*sizeof(long_complex_t));
}


void ffluid_math_get_hamiltonian(data_ptr in, long_double_t *hamiltonian, long_double_t *surf_energy) {
  /* get hamiltonian (constant of motion) */
  unsigned long N = in->N;

  /**/
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[0][j] = -1.IL*(in->V[j]*in->du[j]/in->R[j])/N;
  }
  fftwl_execute(FFTLocal.fp[0]);
  memset(AuxLocal.Y[0]+N/2, 0, N/2);
  *hamiltonian = 0.0L;
  for (long int j = N/2-1; j > -1; j--) {
    *hamiltonian += AuxLocal.Y[0][j]*conjl(AuxLocal.Y[0][j])/( j+1 );
  }
}

void ffluid_math_get_hamiltonian_RV(data_ptr in, long_double_t *hamiltonian, long_double_t *surf_energy) {
}

void ffluid_math_get_angular(data_ptr in, long_double_t *angular_m) {
}
void ffluid_math_get_angular_RV(data_ptr in, long_double_t *angular_m) {
}

void ffluid_math_set_zero_mode(data_ptr in, long_complex_t S0, long_complex_t *out) {
  unsigned long		N = in->N;
  long_complex_t	w = cexpl(1.IL*in->q0);
  long_complex_t	tmp;
  long_double_t 	b = 0.5L*(1.L + powl(in->l, 2))/in->l;
  long_double_t 	xi = (1.L - powl(in->l, 2))/(1.L + powl(in->l, 2));

  AuxLocal.X[0][0] = -0.5L*xi*conjl(w);
  for (unsigned long j = 1; j < N/2 - 1; j++) {
    AuxLocal.X[0][j] = AuxLocal.X[0][0]/(1.L - conjl(AuxLocal.X[0][0])*AuxLocal.X[0][j-1]);
  }
  AuxLocal.X[1][0] = in->R[1]/b - 2.0L*S0*conjl(AuxLocal.X[0][0]);
  for (unsigned long j = 1; j < N/2 - 1; j++) {
    AuxLocal.X[1][j] = in->R[j+1]/b - conjl(AuxLocal.X[0][0])*AuxLocal.X[1][j-1];
    AuxLocal.X[1][j] = AuxLocal.X[1][j]/(1.0L - conjl(AuxLocal.X[0][0])*AuxLocal.X[0][j-1]);
  }

  tmp = AuxLocal.X[1][N/2-2];
  for (long int j = N/2-3; j > -1; j--) {
    tmp = AuxLocal.X[1][j]-AuxLocal.X[0][j]*tmp;
  }
  *out = b*(tmp*AuxLocal.X[0][0] + 1.0L*S0);
}

void ffluid_math_get_surface_spectrum(data_ptr in, data_ptr out) {
  /* get surface potential \Phi(q) and z(q), but the zero modes are not set. */
  /* auxiliary variable S = z_u u_q is introduced */
  unsigned long 	N = in->N;
  
  /* direct inverse O(N log N) */
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[0][j] = in->R[j]/N;
    AuxLocal.X[1][j] = in->V[j]/N;
  }
  /* can be improved with advanced FFTW interface */
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  memcpy(out->R, AuxLocal.Y[0], N*sizeof(long_complex_t));
  memcpy(out->V, AuxLocal.Y[1], N*sizeof(long_complex_t));
}

void ffluid_math_get_volume_RV(data_ptr in, long_double_t *volume) {}

