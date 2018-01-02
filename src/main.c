#include "config.h"
#include "ffluid.h"

control_params 	Control;
evolve_params 	EvolveConfig;
sim_data 	DataCurr, DataPrev;
grid		GridCurr, GridPrev;

/* default grid structure */
const grid GRID_DEFAULT = { false, 0, NULL, NULL };

void ffluid_list_modules() {
  ffluid_io_module();
  ffluid_memory_module();
  ffluid_arrayfunc_module();
  ffluid_mapping_module();
  ffluid_equations_module();
  ffluid_timemarching_module();
}

int main (int argc, char **argv) {
  Control.DataPtrCurr = &DataCurr;
  Control.DataPtrPrev = &DataPrev;
  Control.EvolvePtr = &EvolveConfig;
  GridCurr = GRID_DEFAULT;
  GridPrev = GRID_DEFAULT;

  /* list all functions with their locations */
  ffluid_list_modules();

  ffluid_read_cl_arguments(argc, argv);
  ffluid_set_initial_data();
  ffluid_write_surface(&DataCurr, "test.file");
  //ffluid_evolve();
  return 0;
}
