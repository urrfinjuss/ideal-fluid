#include "ffluid.h"

void init_stepping(double final_time, double max_cfl) {
  double remaining_time = final_time - Data.time;
  double max_dt = max_cfl*pow(2.*PI/(Data.N), 2);

  Evolve_Config.nsteps 	= floor(remaining_time/max_dt) + 1;
  Evolve_Config.dt	= remaining_time/Evolve_Config.nsteps;
  Evolve_Config.cfl	= Evolve_Config.dt*pow(2.*PI/Data.N, -2);
}
