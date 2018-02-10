#include "ffluid.h"

static char surf_name[80];
static char spec_name[80];

void ffluid_timemarching_module() {
  printf("Module ffluid/timemarching.h:\n");
  printf("stepping.c:\tffluid_timemarching_module ffluid_init_stepping ffluid_evolve\n");
  printf("runge_kutta_4.c:\tffluid_runge_kutta_4\n");
  printf("runge_kutta_6.c:\tffluid_runge_kutta_6\n");
}

void ffluid_setup_stepping() {
  __float128 remaining_time = EvolveConfig.final_time - DataCurr.time;
  __float128 max_dt = EvolveConfig.max_cfl*powq(2.Q*PI/(DataCurr.N), 2);

  EvolveConfig.cur_step	= 0;
  EvolveConfig.nsteps 	= floorq(remaining_time/max_dt) + 1;
  EvolveConfig.dt	= remaining_time/EvolveConfig.nsteps;
  EvolveConfig.cfl	= EvolveConfig.dt*powq(2.0Q*M_PIq/DataCurr.N, -2);
  EvolveConfig.dmp_cnt  = 0;
}

void ffluid_detailed_snapshot() {
  if ((EvolveConfig.cur_step % 128) == 0) {
    sprintf(surf_name, "disc/surf_%04lu.txt", EvolveConfig.dmp_cnt);
    sprintf(spec_name, "disc/spec_%04lu.txt", EvolveConfig.dmp_cnt);
    ffluid_math_get_volume(&DataCurr, &DataCurr.Volume);
    ffluid_math_get_surface_variables(&DataCurr, &DataSurface);
    ffluid_math_get_surface_spectrum(&DataCurr, &DataSpectrum);
    printf("Writing data/disc at time %.8Qe r0 = %.15Le\tVolume = %.15Le\n", DataCurr.time, DataCurr.r0, DataCurr.Volume);
    ffluid_write_surface(&DataSurface, surf_name);
    ffluid_write_spectrum(&DataSpectrum, spec_name);
    EvolveConfig.dmp_cnt++; 
  }
}

void ffluid_evolve() {
  ffluid_math_get_volume(&DataCurr, &DataCurr.Volume);
  //ffluid_math_get_r0(&DataCurr);
  ffluid_detailed_snapshot();
  while (EvolveConfig.cur_step < EvolveConfig.nsteps) {
    ffluid_data_copy(&DataCurr, &DataPrev);
    ffluid_runge_kutta_4(&DataCurr);
    //ffluid_math_get_r0(&DataCurr); // this must be part of RK4/6 
    //ffluid_mapping_test_resolved(fname);
    ffluid_detailed_snapshot();
    if (EvolveConfig.cur_step == 32768) break;
  }
  ffluid_detailed_snapshot();
  
  //ffluid_math_get_volume(&DataCurr, &DataCurr.Volume);
}
