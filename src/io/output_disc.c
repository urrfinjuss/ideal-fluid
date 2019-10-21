#include "ffluid.h"


void ffluid_write_surface(data_ptr in, char *fname) {
  unsigned long N = in->N;
  __float128	time = in->time;
  char full_path[160];

  // 
  strcpy(full_path, Control.data_path);
  strcat(full_path, fname);
  //printf("Mapping: u0,q0,l = %.4Le, %.4Le, %4Le\n", in->u0, in->q0, in->l);
  FILE *fh = fopen(full_path, "w");
  fprintf(fh, "# 1. u 2.-3. x+iy  4.-5. psi_u + iHpsi_u\n");
  fprintf(fh, "# Time = %.12Qe\tus = %.16Le\tqs = %.16Le\tl = %.16Le\n\n", DataCurr.time, in->u0, in->q0, in->l);
  for (unsigned long j = 0; j < N; j++) {
    fprintf(fh, "%.16LE\t", in->q[j]);
    fprintf(fh, "%.16LE\t%.16LE\t", creall(in->R[j]), cimagl(in->R[j]));
    fprintf(fh, "%.16LE\t%.16LE\n", creall(in->V[j]), cimagl(in->V[j]));
  }
  fclose(fh);
}

void ffluid_write_full_data(data_ptr inZPh, data_ptr inRV, char *fname) {
  unsigned long N = inZPh->N;
  __float128	time = inZPh->time;
  char full_path[160];

  // 
  strcpy(full_path, Control.data_path);
  strcat(full_path, fname);
  //printf("Mapping: u0,q0,l = %.4Le, %.4Le, %4Le\n", in->u0, in->q0, in->l);
  FILE *fh = fopen(full_path, "w");
  fprintf(fh, "# 1. u 2.-3. R  4.-5. psi_u + iHpsi_u 6.-7. V\n");
  fprintf(fh, "# Time = %.12Qe\tus = %.16Le\tqs = %.16Le\tl = %.16Le\n\n", DataCurr.time, inZPh->u0, inZPh->q0, inZPh->l);
  for (unsigned long j = 0; j < N; j++) {
    fprintf(fh, "%.16LE\t", inZPh->q[j]);
    fprintf(fh, "%.16LE\t%.16LE\t", creall(inRV->R[j]), cimagl(inRV->R[j]));
    fprintf(fh, "%.16LE\t%.16LE\t", creall(inZPh->V[j]), cimagl(inZPh->V[j]));
    fprintf(fh, "%.16LE\t%.16LE\n", creall(inRV->V[j]), cimagl(inRV->V[j]));
  }
  fclose(fh);
}

void ffluid_write_spectrum(data_ptr in, char *fname) {
  unsigned long N = in->N;
  long int k;
  char full_path[160];

  strcpy(full_path, Control.data_path);
  strcat(full_path, fname);
  FILE *fh = fopen(full_path, "w");
  fprintf(fh, "# 1. k 2. |R_k| 4. |V_k|\n");
  fprintf(fh, "# Time = %.12Qe\tus = %.16Le\tqs = %.16Le\tl = %.16Le\n\n", DataCurr.time, in->u0, in->q0, in->l);
  for (unsigned long j = 0; j < N; j++) {
    k = -j;
    if (j > N/2) k += N;
    fprintf(fh, "%24.16LE\t", 1.L*k);
    fprintf(fh, "%24.16LE\t%.16LE\n", cabsl(in->R[j]), cabsl(in->V[j]));
  }
  fclose(fh);
}

void ffluid_write_full_spectrum(data_ptr inRV, data_ptr inZPh, char *fname) {
  unsigned long N = inRV->N;
  long int k;
  char full_path[160];

  strcpy(full_path, Control.data_path);
  strcat(full_path, fname);
  FILE *fh = fopen(full_path, "w");
  fprintf(fh, "# 1. k 2. |R_k| 3. |V_k| 4. |z_k| 5. |F(Phi_u)|\n");
  fprintf(fh, "# Time = %.12Qe\tus = %.16Le\tqs = %.16Le\tl = %.16Le\n\n", DataCurr.time, inRV->u0, inRV->q0, inRV->l);
  for (unsigned long j = 0; j < N; j++) {
    k = -j;
    if (j > N/2) k += N;
    fprintf(fh, "%24.16LE\t", 1.L*k);
    fprintf(fh, "%24.16LE\t%.16LE\t", cabsl(inRV->R[j]), cabsl(inRV->V[j]));
    fprintf(fh, "%24.16LE\t%.16LE\n", cabsl(inZPh->R[j]), cabsl(inZPh->V[j]));
  }
  fclose(fh);
}

void ffluid_start_log(char *fname) {
  char full_path[160];

  strcpy(full_path, Control.data_path);
  strcat(full_path, fname);
  FILE *fh = fopen(full_path, "w");
  fprintf(fh, "# 1. Time 2. y(-pi) 3. Volume (const) 4. KEnergy 5. SEnergy 6. Angular 7. z0\n\n");
  fclose(fh);
}

void ffluid_append_to_log(data_ptr in, char *fname) {
  unsigned long N = in->N;
  char full_path[160];

  strcpy(full_path, Control.data_path);
  strcat(full_path, fname);
  FILE *fh = fopen(full_path, "a");
  //fprintf(fh, "%.12Qe\t%.16Le\t%.16Le\t%.16Le\t%.16Le\n", DataCurr.time, cimagl(in->R[0]), DataCurr.Volume, DataCurr.Hamiltonian, DataCurr.SurfaceEnergy);
  fprintf(fh, "%.12Qe\t%.16Le\t", DataCurr.time, cimagl(in->R[0]));
  fprintf(fh, "%.16Le\t", DataCurr.Volume);
  fprintf(fh, "%.16Le\t", DataCurr.Hamiltonian); 
  fprintf(fh, "%.16Le\t", DataCurr.SurfaceEnergy); 
  fprintf(fh, "%.16Le\t", DataCurr.AngularM); 
  fprintf(fh, "%.16Le\t", cimagl(DataCurr.z0)); 
  fprintf(fh, "%.16Le\n", creall(DataCurr.z0)); 
  fclose(fh);

}
