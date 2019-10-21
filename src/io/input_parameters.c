#include "ffluid.h"

/* Set number of parameters to be set in the input file */
#define INPUT_FILE_ARGS 7	

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
    } else if (!strcmp(keyword,"sft_coeff")) {
      printf("Found:\t%s\n", keyword);
      Control.Sigma = strtoflt128(str_value, NULL);
      counter++;
    } else if (!strcmp(keyword,"fin_time")) {
      printf("Found:\t%s\n", keyword);
      EvolveConfig.final_time = strtoflt128(str_value, NULL);
      counter++;
    } else if (!strcmp(keyword,"number_points")) {
      printf("Found:\t%s\n", keyword);
      DataCurr.N = strtol(str_value, NULL, 10);
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



