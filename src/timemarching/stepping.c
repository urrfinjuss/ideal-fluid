#include "ffluid.h"


void ffluid_timemarching_module() {
  printf("Module ffluid/timemarching.h:\n");
  printf("stepping.c:\tffluid_timemarching_module ffluid_init_stepping ffluid_evolve\n");
  printf("runge_kutta_4.c:\tffluid_runge_kutta_4\n");
  printf("runge_kutta_6.c:\tffluid_runge_kutta_6\n");
}

void ffluid_init_stepping() {
  __float128 remaining_time = EvolveConfig.final_time - DataCurr.time;
  __float128 max_dt = EvolveConfig.max_cfl*powq(2.Q*PI/(DataCurr.N), 2);

  EvolveConfig.cur_step	= 0;
  EvolveConfig.nsteps 	= floorq(remaining_time/max_dt) + 1;
  EvolveConfig.dt	= remaining_time/EvolveConfig.nsteps;
  EvolveConfig.cfl	= EvolveConfig.dt*powq(2.0Q*M_PIq/DataCurr.N, -2);

  ffluid_init_runge_kutta_4();
  ffluid_data_init_copy(&DataCurr, &DataPrev);
}

void ffluid_evolve() {
  char 		fname[80];
  long_double_t	Volume;
  printf("\tffluid_init_stepping()\n");
  ffluid_init_stepping();
  while (EvolveConfig.cur_step < EvolveConfig.nsteps) {
    printf("\titer = %lu\n", EvolveConfig.cur_step);
    ffluid_data_copy(&DataCurr, &DataPrev);
    ffluid_runge_kutta_4(&DataCurr);
    ffluid_mapping_test_resolved(fname);

    sprintf(fname, "disc_surf_%04lu.txt", EvolveConfig.cur_step);
    ffluid_math_get_surface_variables(&DataCurr, &DataSurface);
    ffluid_write_surface(&DataSurface, fname);
    //ffluid_math_get_volume(&DataCurr, &Volume);

    if (EvolveConfig.cur_step == 10) break;
    //ffluid_write_surface();
    //ffluid_append_to_log();
  }
}
