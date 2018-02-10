#include "config.h"
#include "ffluid.h"

control_params 	Control;
evolve_params 	EvolveConfig;
sim_data 	DataCurr, DataPrev;
sim_data 	DataSpectrum, DataSurface;

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

  /* ffluid_topread() */
  ffluid_read_cl_arguments(argc, argv);
  ffluid_set_initial_data(&DataCurr);

  /* ffluid_alloc_timemarching() */
  ffluid_data_init_copy(&DataCurr, &DataPrev);
  ffluid_init_runge_kutta_4(); // or rk6
  ffluid_alloc_equations();

  /* ffluid_alloc_output() */
  ffluid_data_init_copy(&DataCurr, &DataSpectrum); 
  ffluid_data_init_copy(&DataCurr, &DataSurface); 
  ffluid_math_init_surface();

  /* ffluid_timemarching() */
  ffluid_setup_stepping();
  ffluid_evolve();
  printf("Complete\n");
  /* unused block for output */
  //ffluid_data_init_copy(&DataCurr, &DataSurface);
  //ffluid_math_init_surface();
  //ffluid_math_get_surface_variables(&DataCurr, &DataSurface);
  //ffluid_write_surface(&DataSurface, "test.disc.file");
  //ffluid_math_get_volume(&DataCurr, &Volume);

  return 0;
}
