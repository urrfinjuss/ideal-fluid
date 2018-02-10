#include "ffluid.h"

void ffluid_set_initial_data(data_ptr in) {
  /* either read initial data from file, or generate new initial data? */
  ffluid_init_data(in);
  ffluid_read_initial_data(in);
  ffluid_setup_grid(in);
}

void ffluid_read_initial_data(data_ptr in) {
  FILE *fh = fopen(Control.res_name,"r");
  char line[512];
  
  if (fh) {
    ffluid_read_mapping_parameters(in, fh);
    ffluid_read_data_from_file(in, fh);
  } else {
    printf("File %s not found.\n", Control.res_name);
    exit(0);
  }
}

void ffluid_read_mapping_parameters(data_ptr in, FILE *fh) {
  char line[512], *val[4];
  for (int j = 0; j < 4; j++) val[j] = malloc(80*sizeof(char));
  if (fgets(line, 512, fh));
  if (fgets(line, 512, fh));
  sscanf(line, "# Time = %s\tus = %s\tqs = %s\tl = %s", val[0], val[1], val[2], val[3]); 
  in->time = strtoflt128(val[0], NULL);
  in->u0 = strtold(val[1], NULL);
  in->q0 = strtold(val[2], NULL);
  in->l = strtold(val[3], NULL);
  for (int j = 3; j > -1; j--) free(val[j]);
}

void ffluid_read_data_from_file(data_ptr in, FILE *fh) {
  unsigned long counter = 0;
  char line[512];
  char *val[5];

  if (fgets(line, 512, fh));
  for (int j = 0; j < 5; j++) val[j] = malloc(80*sizeof(char));
  while ((counter < DataCurr.N)&&(fgets(line, 512, fh))) {
    sscanf(line, "%s\t%s\t%s\t%s\t%s\n", val[0], val[1], val[2], val[3], val[4]);
    in->R[counter] = strtold(val[1], NULL) + 1.0IL*strtod(val[2], NULL);
    in->V[counter] = strtold(val[3], NULL) + 1.0IL*strtod(val[4], NULL);
    counter++;
  }
  for (int j = 4; j > -1; j--) free(val[j]);
}


