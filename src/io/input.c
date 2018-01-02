#include "ffluid.h"

/* Set number of parameters to be set in the input file */
#define INPUT_FILE_ARGS 5	

void ffluid_scan_input_file(FILE *fh) {
  unsigned int nparameters = INPUT_FILE_ARGS;
  unsigned int counter = 0;
  char line[80], keyword[80], str_value[80];

  while ((counter < nparameters)&&(fgets(line, 80, fh))) {
    sscanf(line, "%s\t%s", keyword, str_value);
    if (!strcmp(keyword,"run_name")) {
      printf("Found:\t%s\n", keyword);
      sprintf(Control.run_name, "%s", str_value);
      counter++;
    } else if (!strcmp(keyword,"res_name")) {
      printf("Found:\t%s\n", keyword);
      sprintf(Control.res_name, "%s", str_value);
      counter++;
    } else if (!strcmp(keyword,"data_path")) {
      printf("Found:\t%s\n", keyword);
      sprintf(Control.data_path, "%s", str_value);
      counter++;
    } else if (!strcmp(keyword,"cfl_value")) {
      printf("Found:\t%s\n", keyword);
      EvolveConfig.max_cfl = strtoflt128(str_value, NULL);
      counter++;
    } else if (!strcmp(keyword,"number_points")) {
      printf("Found:\t%s\n", keyword);
      DataCurr.N = strtol(str_value, NULL, 10);
      DataPrev.N = strtol(str_value, NULL, 10);
      counter++;
    } 
  }
  if (counter == nparameters) {
    printf("Scan Complete: %d lines\n", counter);
  } else {
    printf("Configuration file incomplete.\n");
    exit(0);
  }
}

void ffluid_read_input_file(char *filename) {
  FILE *fh = fopen(filename,"r");
  if (fh) {
    printf("File %s opened.\nVerifying set parameters:\n", filename);
    ffluid_scan_input_file(fh);
    fclose(fh);
  } else { 
    printf("File %s not found.\n", filename);
    exit(0);
  }
}

void ffluid_read_cl_arguments(int narg, char **argv) {
  if (narg != 2) {
    display_help();
    exit(0);
  } else {
    if (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-h")) {
      display_help();
      exit(0);
    } else if (!strcmp(argv[1], "--version") || !strcmp(argv[1],"-v")) {
      display_version(); 
      exit(0);
    } else {
      ffluid_read_input_file(argv[1]);
    }
  }
}

void ffluid_read_mapping_parameters(FILE *fh) {
  char line[512], *val[4];
  for (int j = 0; j < 4; j++) val[j] = malloc(80*sizeof(char));
  if (fgets(line, 512, fh));
  if (fgets(line, 512, fh));
  sscanf(line, "# Time = %s\tus = %s\tqs = %s\tl = %s", val[0], val[1], val[2], val[3]); 
  DataCurr.time = strtoflt128(val[0], NULL);
  DataCurr.u0 = strtod(val[1], NULL);
  DataCurr.q0 = strtod(val[2], NULL);
  DataCurr.l = strtod(val[3], NULL);
  for (int j = 3; j > -1; j--) free(val[j]);
}

void ffluid_read_data_from_file(FILE *fh) {
  unsigned long counter = 0;
  char line[512];
  char *val[5];

  if (fgets(line, 512, fh));
  for (int j = 0; j < 5; j++) val[j] = malloc(80*sizeof(char));
  while ((counter < DataCurr.N)&&(fgets(line, 512, fh))) {
    sscanf(line, "%s\t%s\t%s\t%s\t%s\n", val[0], val[1], val[2], val[3], val[4]);
    DataCurr.Q[counter] = strtod(val[1], NULL) + 1.0I*strtod(val[2], NULL);
    DataCurr.V[counter] = strtod(val[3], NULL) + 1.0I*strtod(val[4], NULL);
    counter++;
  }
  for (int j = 4; j > -1; j--) free(val[j]);
}

void ffluid_read_initial_data() {
  FILE *fh = fopen(Control.res_name,"r");
  char line[512];
  
  if (fh) {
    ffluid_read_mapping_parameters(fh);
    ffluid_read_data_from_file(fh);
  } else {
    printf("File %s not found.\n", Control.res_name);
    exit(0);
  }
}

void ffluid_set_initial_data() {
  /* either read initial data from file, or generate new initial data? */
  ffluid_init_data(&DataCurr, DataCurr.N);
  ffluid_read_initial_data();
}

