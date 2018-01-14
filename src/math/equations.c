#include "ffluid.h"

static sim_data SimLocal;
static aux_data AuxLocal;
static fft_list FFTLocal;

void ffluid_math_init_equations() {
  ffluid_data_init_copy(&DataCurr, &SimLocal);
  ffluid_alloc_aux_array(&AuxLocal, 5, DataCurr.N);
  ffluid_alloc_fft_plans(&AuxLocal, &FFTLocal);
}

void ffluid_equations_module() {
  printf("Module ffluid/equations.h:\n");
  printf("equations.c:\tffluid_arrayfunc_module ffluid_call_rhs\n");
}

/* 
void ffluid_projection(data_ptr in) {
  unsigned long N = in->N;
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[0][j] = 2.0L*creall(in->Q[j]*in->Q[j]*conjl(in->V[j]))/N;
    AuxLocal.X[1][j] = creall(in->V[j]*conjl(in->V[j]))/N;
  }
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  memset(AuxLocal.Y[0]+N/2, 0, N/2);
  memset(AuxLocal.Y[1]+N/2, 0, N/2);
  AuxLocal.Y[0][0] = 0.5L*AuxLocal.Y[0][0];
  AuxLocal.Y[1][0] = 0.5L*AuxLocal.Y[1][0];
  fftwl_execute(FFTLocal.bp[0]);
  fftwl_execute(FFTLocal.bp[1]);
  memcpy(SimLocal.Q, AuxLocal.Y[0], N*sizeof(long_complex_t));
  memcpy(SimLocal.V, AuxLocal.Y[1], N*sizeof(long_complex_t));
}
*/

void ffluid_call_rhs(data_ptr in, data_ptr out) {
  unsigned long 	N = in->N;
  long_complex_t	w1 = cexpl( 1.IL*in->u0 - 2.0L*atanhl(in->l));
  long_complex_t	w2 = cexpl(-1.IL*in->u0 - 2.0L*atanhl(in->l));
  long_complex_t	b1U = 0.0L, b2U = 0.0L;

  if (in->l == 1.0L) {
    w2 = 0.0L;
    w1 = 0.0L;
  }
  memcpy(AuxLocal.X[0], in->Q, N*sizeof(long_complex_t));
  memcpy(AuxLocal.X[1], in->V, N*sizeof(long_complex_t));
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  memset(AuxLocal.Y[0]+N/2, 0, N/2);
  memset(AuxLocal.Y[1]+N/2, 0, N/2);
  for (unsigned long j = 0; j < N/2; j++) {
    AuxLocal.Y[0][j] = -1.IL*j*AuxLocal.Y[0][j]/N;
    AuxLocal.Y[1][j] = -1.IL*j*AuxLocal.Y[1][j]/N;
  }
  fftwl_execute(FFTLocal.bp[0]);
  fftwl_execute(FFTLocal.bp[1]);
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[2][j]= 2.0L*creall(in->V[j]*conjl(in->Q[j]*in->Q[j]))/N;
    AuxLocal.X[3][j]= in->V[j]*conjl(in->V[j])/N;
  }
  fftwl_execute(FFTLocal.fp[2]);
  fftwl_execute(FFTLocal.fp[3]);
  AuxLocal.Y[4][0] = 0.0L;
  b1U = AuxLocal.Y[2][N/2+1];
  b2U = AuxLocal.Y[2][N/2-1];
  for (unsigned long j = 1; j < N/2 - 1; j++) {
    b1U = b1U*w1 + AuxLocal.Y[2][N/2+1+j];
    b2U = b2U*w2 + AuxLocal.Y[2][N/2-1-j];
  }
  memset(AuxLocal.Y[2]+N/2, 0, N/2*sizeof(long_complex_t));
  memset(AuxLocal.Y[3]+N/2, 0, N/2*sizeof(long_complex_t));
  memset(AuxLocal.Y[4]+N/2, 0, N/2*sizeof(long_complex_t));
  AuxLocal.Y[2][0] = 0.5L*AuxLocal.Y[2][0];
  for (unsigned long j = 0; j < N/2; j++) {
    AuxLocal.Y[3][j] = -1.0IL*j*AuxLocal.Y[3][j];
    AuxLocal.Y[4][j] = -1.0IL*j*AuxLocal.Y[2][j];
  }
  fftwl_execute(FFTLocal.bp[2]);
  fftwl_execute(FFTLocal.bp[3]);
  fftwl_execute(FFTLocal.bp[4]);
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[2][j] += 0.5L*(b1U*w1 - b2U*w2);
  }
  for (long int j = 0; j < N; j++) {
    out->Q[j] = 0.5IL*(2.L*AuxLocal.X[0][j]*AuxLocal.X[2][j]-AuxLocal.X[4][j]*in->Q[j]);
    out->V[j] = 1.0IL*(AuxLocal.X[2][j]*AuxLocal.X[1][j]-in->Q[j]*in->Q[j]*AuxLocal.X[3][j]);
    out->Q[j] = out->Q[j]*in->du[j];
    out->V[j] = out->V[j]*in->du[j];
  }
}





























