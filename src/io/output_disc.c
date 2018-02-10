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
  fprintf(fh, "# 1. u 2. Q.re, Q.im 4. V.re, V.im\n");
  fprintf(fh, "# Time = %Qe\tus = %.16Le\tqs = %.16Le\tl = %.16Le\n\n", in->time, in->u0, in->q0, in->l);
  for (unsigned long j = 0; j < N; j++) {
    fprintf(fh, "%.16LE\t", in->q[j]);
    fprintf(fh, "%.16LE\t%.16LE\t", creall(in->R[j]), cimagl(in->R[j]));
    fprintf(fh, "%.16LE\t%.16LE\n", creall(in->V[j]), cimagl(in->V[j]));
  }
  fclose(fh);
}

void ffluid_write_spectrum(data_ptr in, char *fname) {
  unsigned long N = in->N;
  char full_path[160];

  strcpy(full_path, Control.data_path);
  strcat(full_path, fname);
  FILE *fh = fopen(full_path, "w");
  fprintf(fh, "# 1. k 2. |Q_k| 4. |V_k|\n");
  fprintf(fh, "# Time = %Qe\tus = %.16Le\tqs = %.16Le\tl = %.16Le\n\n", in->time, in->u0, in->q0, in->l);
  for (unsigned long j = 0; j < N/2; j++) {
    fprintf(fh, "%.16LE\t", -1.0L*j);
    fprintf(fh, "%.16LE\t%.16LE\n", cabsl(in->R[j]), cabsl(in->V[j]));
  }
  fclose(fh);
}

void ffluid_append_to_log() {
}
