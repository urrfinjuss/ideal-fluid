#include "ffluid.h"


void ffluid_timemarching_module() {
  printf("Module ffluid/timemarching.h:\n");
  printf("stepping.c:\tffluid_timemarching_module ffluid_init_stepping ffluid_evolve\n");
  printf("runge_kutta_4.c:\tffluid_runge_kutta_4\n");
  printf("runge_kutta_6.c:\tffluid_runge_kutta_6\n");
}

void ffluid_init_stepping() {
  double remaining_time = EvolveConfig.final_time - DataCurr.time;
  double max_dt = EvolveConfig.max_cfl*pow(2.Q*PI/(DataCurr.N), 2);

  EvolveConfig.cur_step	= 0;
  EvolveConfig.nsteps 	= floorq(remaining_time/max_dt) + 1;
  EvolveConfig.dt	= remaining_time/EvolveConfig.nsteps;
  EvolveConfig.cfl	= EvolveConfig.dt*powq(2.Q*M_PIq/DataCurr.N, -2);
}

void ffluid_evolve() {
  ffluid_init_stepping();
  while (EvolveConfig.cur_step < EvolveConfig.nsteps) {
    ffluid_data_copy(&DataCurr, &DataPrev);
    ffluid_runge_kutta_4(&DataCurr);
    ffluid_mapping_test_resolved();
    //ffluid_write_surface();
    //ffluid_append_to_log();
  }
}
