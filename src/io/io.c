#include "ffluid.h"

void ffluid_io_module() {
  printf("Module ffluid/io.h:\n");
  printf("io.c:\t\tffluid_io_module\n");
  printf("input.c:\tffluid_scan_input_file ffluid_read_input_file ffluid_read_cl_arguments ffluid_read_initial_data fflui_set_initial_data ffluid_read_mapping_parameters ffluid_read_data_from_file\n");
  printf("output.c:\tffluid_write_surface ffluid_append_to_log\n");
}

void ffluid_write_array(long_complex_t *in, unsigned long N, char *fname) {
  char full_path[160];

  strcpy(full_path, Control.data_path);
  strcat(full_path, fname);
  FILE *fh = fopen(full_path, "w");
  fprintf(fh, "# 1. u 2. X 3. Y.\n\n");
  for (unsigned long j = 0; j < N; j++) {
    fprintf(fh, "%.16LE\t", PI*(2.0L*j/N - 1.L));
    fprintf(fh, "%.16LE\t%.16LE\n", creall(in[j]), cimagl(in[j]));
  }
  fclose(fh);
}
