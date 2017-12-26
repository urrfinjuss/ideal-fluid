#include "ffluid.h"

void test_input_file(char *filename) {
  FILE *fh = fopen(filename,"r");
  if (fh) {
    printf("File %s opened.\nChecking set parameters:\nSTUB", filename);
  } else { 
    printf("File %s not found.\n", filename);
    exit(0);
  }
}

void ffluid_read_cl_arguments(int narg, char **argv) {
  printf("narg = %d\n", narg);
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
      test_input_file(argv[1]);
    }
  }
}

