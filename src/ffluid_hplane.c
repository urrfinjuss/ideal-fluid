#include "config.h"
#include "ffluid.h"

control_params 	Control;
evolve_params 	EvolveConfig;
sim_data 	DataCurr, DataPrev, DataSurface;

/* list functions */
void ffluid_list_modules() {
  ffluid_io_module();
  ffluid_memory_module();
  ffluid_arrayfunc_module();
  ffluid_mapping_module();
  ffluid_equations_module();
  ffluid_timemarching_module();
}

/* main function */
int main (int argc, char **argv) {
  Control.DataPtrCurr = &DataCurr;
  Control.DataPtrPrev = &DataPrev;
  Control.EvolvePtr = &EvolveConfig;

  /* list all functions with their locations */
  // ffluid_list_modules();

  ffluid_read_cl_arguments(argc, argv);
  ffluid_set_initial_data(&DataCurr);
  
  ffluid_data_init_copy(&DataCurr, &DataSurface);
  printf("ffluid_data_init_copy()\n");

  ffluid_math_init_surface();
  printf("ffluid_init_surface_math()\n");
  ffluid_math_get_surface_variables(&DataCurr, &DataSurface);
  printf("ffluid_get_surface_variables()\n");
  ffluid_write_surface(&DataSurface, "test.hplane.file");
  
  //ffluid_evolve();
  return 0;
}
