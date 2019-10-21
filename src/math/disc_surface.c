#include "ffluid.h"

static sim_data SimLocal;
static aux_data AuxLocal;
static fft_list FFTLocal;

void ffluid_math_init_surface() {
  ffluid_data_init_copy(&DataCurr, &SimLocal);
  ffluid_alloc_aux_array(&AuxLocal, 2, DataCurr.N);
  ffluid_alloc_fft_plans(&AuxLocal, &FFTLocal);
}

void ffluid_math_get_volume(data_ptr in, long_double_t *volume){
  /* get fluid volume (constant of motion) */
  unsigned long N = in->N;

  /* direct inverse O(N log N) */
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[0][j] = (cpowl(in->R[j], -2))/N;
  }
  fftwl_execute(FFTLocal.fp[0]);
  memset(AuxLocal.Y[0]+N/2, 0, N/2);
  *volume = 0.0L;
  for (long int j = N/2-1; j > -1; j--) {
    *volume += AuxLocal.Y[0][j]*conjl(AuxLocal.Y[0][j])/(j + 1);
  }
  *volume = PI*(*volume);
}

void ffluid_math_get_hamiltonian(data_ptr in, long_double_t *hamiltonian, long_double_t *surf_energy) {
  /* get hamiltonian (constant of motion) */
  unsigned long N = in->N;
  /**/
  for (unsigned long j = 0; j < N; j++) {
    //AuxLocal.X[0][j] = -1.IL*(in->V[j]*in->du[j]/in->R[j])/N;
    AuxLocal.X[0][j] = in->V[j]*cpowl(in->R[j], -2)/N;
    AuxLocal.X[1][j] = cpowl(in->R[j], -1)/N;
  }
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  memset(AuxLocal.Y[0]+N/2, 0, N/2);
  *hamiltonian = 0.0L;
  *surf_energy = 4.0L*Control.Sigma*AuxLocal.Y[1][0]*conjl(AuxLocal.Y[1][0]);
  // old way
  /*
  for (long int j = N/2-1; j > -1; j--) {
    *hamiltonian += AuxLocal.Y[0][j]*conjl(AuxLocal.Y[0][j])/( j+1 );
    *surf_energy += 4.0L*Control.Sigma*AuxLocal.Y[1][j]*conjl(AuxLocal.Y[1][j]);
  }
  */
  // by the text
  for (long int j = N/2-1; j > 0; j--) {
    *hamiltonian += 1.0L*AuxLocal.Y[0][j]*conjl(AuxLocal.Y[0][j])/( j );
    *surf_energy += 4.0L*Control.Sigma*AuxLocal.Y[1][j]*conjl(AuxLocal.Y[1][j]);
  }
}

/*
void ffluid_math_get_r0(data_ptr in) {
  unsigned long N = in->N;
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[0][j] = (cpowl(in->R[j], -1))*in->du[j]/N;
  }
  fftwl_execute(FFTLocal.fp[0]);
  
  DataCurr.r0 = 0.L;
  for (long int j = N/2-1; j > 1; j--) {
    DataCurr.r0 += creall(j*AuxLocal.Y[0][j]*conjl(AuxLocal.Y[0][j]));
  }
  //printf("r0 = %.12LE\tVolume = %.12LE\n", DataCurr.r0, DataCurr.Volume);
  DataCurr.r0 = 1.L/sqrtl(DataCurr.Volume/PI - DataCurr.r0);
  //exit(1);
}
*/

/*
void ffluid_math_get_surface_variables(data_ptr in, data_ptr out) { 
  // get surface potential \Phi(q) and z(q), but the zero modes are not set. //
  // auxiliary variable S = z_u u_q is introduced //
  unsigned long 	N = in->N;
  
  // direct inverse O(N log N) //
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[0][j] = (cpowl(in->R[j], -1))*in->du[j]/N;
    AuxLocal.X[1][j] = -1.0IL*in->V[j]*AuxLocal.X[0][j]/N;
  }
  // can be improved with advanced FFTW interface //
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  memset(AuxLocal.Y[0]+N/2, 0, N/2);
  memset(AuxLocal.Y[1]+N/2, 0, N/2);
  for (long int j = N/2; j > -1; j--) {
    // AuxLocal.ft[0] stores Fourier modes of z_k, k < 0   //
    // AuxLocal.ft[1] stores Fourier modes of Phi_k, k < 0 // 
    AuxLocal.Y[0][j+1] = -1.0L*AuxLocal.X[0][j]/(j+1);
    AuxLocal.Y[1][j+1] = -1.0L*AuxLocal.X[1][j]/(j+1);
    //S0 += j*creall(AuxLocal.ft[0][j]*conjl(AuxLocal.ft[0][j]));
    //T0 += -creall(AuxLocal.ft[0][j]);
  }
  AuxLocal.Y[0][0] = 0.0L;
  AuxLocal.Y[1][0] = 0.0L;
  //memcpy(SimLocal.Q, AuxLocal.X[0], N*sizeof(long_complex_t));
  //memcpy(SimLocal.V, AuxLocal.X[1], N*sizeof(long_complex_t));
  //ffluid_write_surface(&SimLocal, "SimLocal2.disc.file");
  // setting of the zero modes of Phi,Z goes here //
  // can be improved with advanced FFTW interface // 
  fftwl_execute(FFTLocal.bp[0]);
  fftwl_execute(FFTLocal.bp[1]);
  memcpy(out->R, AuxLocal.X[0], N*sizeof(long_complex_t));
  memcpy(out->V, AuxLocal.X[1], N*sizeof(long_complex_t));
}
*/

void ffluid_math_get_angular(data_ptr in, long_double_t *angular_m) {
  /* input: surface variables */   
  unsigned long		N = in->N;
  
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[0][j] = (in->R[j])*conjl(in->R[j])/N;
    AuxLocal.X[1][j] = creall(in->V[j])/N;
    //AuxLocal.X[0][j] = in->R[j]/N;
  }
  //ffluid_write_surface(in, "test_potential.txt");
  //exit(0);
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  long_complex_t	AM = 0.L;
  for (unsigned long j = N/2; j > 0; j--) {
    AM +=    -1.0IL*j*AuxLocal.Y[0][j]*conjl(AuxLocal.Y[1][j]);   
    //*AuxLocal.Y[0][j];
    AM +=    -1.0IL*j*AuxLocal.Y[0][N-j]*conjl(AuxLocal.Y[1][N-j]); 
    //*AuxLocal.Y[0][N-j];
  }

  //printf("AM = %.12Le\t%.12Le\n", creall(AM), cimagl(AM)); exit(0);
  *angular_m = -2.L*PI*cimagl(AM);
}

void ffluid_math_get_surface_variables(data_ptr in, data_ptr out) { 
  unsigned long 	N = in->N;
  
  for (unsigned long j = 0; j < N; j++) {
    /* this is dZ in Q */
    AuxLocal.X[0][j] = (cpowl(in->R[j], -2))*in->du[j]/N;
    //AuxLocal.X[1][j] = -1.0IL*in->V[j]*AuxLocal.X[0][j];
    AuxLocal.X[1][j] = -1.0IL*in->V[j]*(cpowl(in->R[j], -2))/N;
  }
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  memset(AuxLocal.Y[0]+N/2, 0, N/2);
  memset(AuxLocal.Y[1]+N/2, 0, N/2);
  for (unsigned long j = N/2; j > 0; j--) {
    AuxLocal.Y[0][j] = -1.0L*AuxLocal.Y[0][j]/j;
    AuxLocal.Y[1][j] = -1.0L*AuxLocal.Y[1][j]/j;
  }
  AuxLocal.Y[0][0] = 0.0L;
  AuxLocal.Y[1][0] = 0.0L;
  fftwl_execute(FFTLocal.bp[0]);
  fftwl_execute(FFTLocal.bp[1]);
  memcpy(out->R, AuxLocal.X[0], N*sizeof(long_complex_t));
  memcpy(out->V, AuxLocal.X[1], N*sizeof(long_complex_t));
  /* we must now set the zero mode of z */
  
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[0][j] = out->R[j]/N;
    AuxLocal.X[1][j] = (cpowl(out->R[j], 2))/N;
  }
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  long_complex_t 	volume = 0.L;
  long_complex_t 	sumsqr = 0.L;
  for (unsigned long j = N/2; j > 0; j--) {
    volume +=  1.0L*j*AuxLocal.Y[0][j]*conjl(AuxLocal.Y[0][j]);
    sumsqr +=  -0.5IL*j*AuxLocal.Y[1][j]*conjl(AuxLocal.Y[0][j]);
  }
  //long_complex_t Z0 = -1.IL*sumsqr/volume;
  AuxLocal.Y[0][0] = -1.IL*sumsqr/volume;
  fftwl_execute(FFTLocal.bp[0]);
  memcpy(out->R, AuxLocal.X[0], N*sizeof(long_complex_t));
  //printf("Vol = %.12Le\n", volume);
  //printf("SmS = %.12Le\t%.12Le\n", creall(Z0), cimagl(Z0));
  //exit(0); 
}



// ---- Following the notation of the text main_v06.pdf  ---- //

void ffluid_math_get_surface_spectrum(data_ptr in, data_ptr out) {
  /* get surface potential \Phi_u(q) and z(q), but the zero modes are not set. */
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


void ffluid_math_get_volume_RV(data_ptr in, long_double_t *volume){
  /* get fluid volume (constant of motion) */
  unsigned long N = in->N;

  /* direct inverse O(N log N) */
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[0][j] = (cpowl(in->R[j], -1))/N;
  }
  fftwl_execute(FFTLocal.fp[0]);
  //memset(AuxLocal.Y[0]+N/2, 0, N/2);
  *volume = cabsl(AuxLocal.Y[0][N/2-1]);
  //*volume = 0.0L;
  for (long int j = N/2-1; j > 0; j--) {
    *volume += AuxLocal.Y[0][j]*conjl(AuxLocal.Y[0][j])/j;
  }
  *volume = PI*(*volume);
}

void ffluid_math_get_hamiltonian_RV(data_ptr in, long_double_t *hamiltonian, long_double_t *surf_energy) {
  /* get hamiltonian (constant of motion) */
  unsigned long N = in->N;
  /**/
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[0][j] = in->V[j]*cpowl(in->R[j], -1)/N;
    AuxLocal.X[1][j] =1.L/csqrt(in->du[j]*in->R[j])/N;
  }
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  memset(AuxLocal.Y[0]+N/2, 0, N/2);
  *hamiltonian = 0.0L;
  *surf_energy = AuxLocal.Y[1][0]*conjl(AuxLocal.Y[1][0]);
  for (long int j = N/2-1; j > 0; j--) {
    *hamiltonian += 1.0L*AuxLocal.Y[0][j]*conjl(AuxLocal.Y[0][j])/j;
    *surf_energy += AuxLocal.Y[1][j]*conjl(AuxLocal.Y[1][j]);
  }
  *surf_energy = 4.0L*Control.Sigma*(*surf_energy);
}

void ffluid_math_get_angular_RV(data_ptr in, long_double_t *angular_m) {
  /* input: surface variables */   
  unsigned long		N = in->N;
  
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[0][j] = cabsl(in->R[j])/N;
    AuxLocal.X[1][j] = creall(in->V[j])/N;
  }
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  long_complex_t	AM = 0.L;
  for (unsigned long j = N/2; j > 0; j--) {
    AM +=    -1.0IL*j*AuxLocal.Y[0][j]*conjl(AuxLocal.Y[1][j]);   
    AM +=    -1.0IL*j*AuxLocal.Y[0][N-j]*conjl(AuxLocal.Y[1][N-j]); 
  }
  *angular_m = -2.L*PI*cimagl(AM);
}

void ffluid_math_get_surface_variables_RV(data_ptr in, data_ptr out) { 
  unsigned long 	N = in->N;
  
  for (unsigned long j = 0; j < N; j++) {
    /* this is dZ in R */
    AuxLocal.X[0][j] = (cpowl(in->R[j], -1))/N;
    //AuxLocal.X[1][j] = -1.0IL*in->V[j]*AuxLocal.X[0][j];
    AuxLocal.X[1][j] = -1.0IL*in->V[j]*(cpowl(in->R[j], -1))/N;
  }
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  memset(AuxLocal.Y[0]+N/2, 0, N/2);
  memset(AuxLocal.Y[1]+N/2, 0, N/2);
  for (unsigned long j = N/2; j > 0; j--) {
    AuxLocal.Y[0][j] = 1.0IL*AuxLocal.Y[0][j]/j;  // must be 1.0IL instead of -1.0L
    //AuxLocal.Y[1][j] = -1.0L*AuxLocal.Y[1][j]/j; // we do not integrate \Phi_u for angular momentum
  }
  AuxLocal.Y[0][0] = 0.0L;
  //AuxLocal.Y[1][0] = 0.0L;
  fftwl_execute(FFTLocal.bp[0]);
  fftwl_execute(FFTLocal.bp[1]);
  memcpy(out->R, AuxLocal.X[0], N*sizeof(long_complex_t));
  memcpy(out->V, AuxLocal.X[1], N*sizeof(long_complex_t));
  /* we must now set the zero mode of z */
  
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[0][j] = out->R[j]/N;
    AuxLocal.X[1][j] = (cpowl(out->R[j], 2))/N;
  }
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  long_complex_t 	volume = 0.L;
  long_complex_t 	sumsqr = 0.L;
  for (unsigned long j = N/2; j > 0; j--) {
    volume +=   1.0L*j*AuxLocal.Y[0][j]*conjl(AuxLocal.Y[0][j]);
    sumsqr += -0.5IL*j*AuxLocal.Y[1][j]*conjl(AuxLocal.Y[0][j]);
  }
  //long_complex_t Z0 = -1.IL*sumsqr/volume;
  AuxLocal.Y[0][0] = -1.IL*sumsqr/volume;
  in->z0 = AuxLocal.Y[0][0];
  fftwl_execute(FFTLocal.bp[0]);
  memcpy(out->R, AuxLocal.X[0], N*sizeof(long_complex_t));
  //printf("Vol = %.12Le\n", volume);
  //printf("SmS = %.12Le\t%.12Le\n", creall(Z0), cimagl(Z0));
  //exit(0); 
}
