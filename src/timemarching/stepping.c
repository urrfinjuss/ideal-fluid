#include "ffluid.h"

static char surf_name[80];
static char full_name[80];
static char spec_name[80];

void ffluid_timemarching_module() {
  printf("Module ffluid/timemarching.h:\n");
  printf("stepping.c:\tffluid_timemarching_module ffluid_init_stepping ffluid_evolve\n");
  printf("runge_kutta_4.c:\tffluid_runge_kutta_4\n");
  printf("runge_kutta_6.c:\tffluid_runge_kutta_6\n");
}

void ffluid_setup_stepping() {
  __float128 remaining_time = EvolveConfig.final_time - DataCurr.time;
  __float128 max_dt = EvolveConfig.max_cfl*powq(2.Q/(DataCurr.N), 2);

  EvolveConfig.cur_step	= 0;
  EvolveConfig.nsteps 	= floorq(remaining_time/max_dt) + 1;
  EvolveConfig.dt	= remaining_time/EvolveConfig.nsteps;
  EvolveConfig.cfl	= EvolveConfig.dt*powq(2.0Q/DataCurr.N, -2);
  EvolveConfig.dmp_cnt  = 0;
}
void ffluid_last_detailed_snapshot() {
    //ffluid_math_get_surface_variables(&DataCurr, &DataSurface);
    sprintf(surf_name, "disc/surf_%04lu.txt", EvolveConfig.dmp_cnt);
    sprintf(spec_name, "disc/spec_%04lu.txt", EvolveConfig.dmp_cnt);
    ffluid_math_get_surface_variables(&DataCurr, &DataSurface);
    ffluid_math_get_surface_spectrum(&DataCurr, &DataSpectrum);
    ffluid_math_get_volume(&DataCurr, &DataCurr.Volume);
    ffluid_math_get_hamiltonian(&DataCurr, &DataCurr.Hamiltonian, &DataCurr.SurfaceEnergy);
    ffluid_math_get_angular(&DataSurface, &DataCurr.AngularM);
    /* output with angular momentum */
    ffluid_append_to_log(&DataSurface, "time.log");
    printf("Writing %s at time %.8Qe\t", surf_name, DataCurr.time);
    printf("Vol = %.15Le\t", DataCurr.Volume);
    printf("K+P = %.15Le\t", DataCurr.Hamiltonian+DataCurr.SurfaceEnergy);
    printf("J = %.15Le\n", DataCurr.AngularM);
    ffluid_write_surface(&DataSurface, surf_name);
    ffluid_write_spectrum(&DataSpectrum, spec_name);
    EvolveConfig.dmp_cnt++; 
}

void ffluid_detailed_snapshot() {
  if ((EvolveConfig.cur_step % 64) == 0) {
    ffluid_math_get_surface_variables(&DataCurr, &DataSurface);
    ffluid_math_get_volume(&DataCurr, &DataCurr.Volume);
    ffluid_math_get_hamiltonian(&DataCurr, &DataCurr.Hamiltonian, &DataCurr.SurfaceEnergy);
    ffluid_math_get_angular(&DataSurface, &DataCurr.AngularM);
    ffluid_append_to_log(&DataSurface, "time.log");
  }
  if ((EvolveConfig.cur_step % 128) == 0) {
    sprintf(surf_name, "disc/surf_%04lu.txt", EvolveConfig.dmp_cnt);
    sprintf(spec_name, "disc/spec_%04lu.txt", EvolveConfig.dmp_cnt);
    ffluid_math_get_surface_variables(&DataCurr, &DataSurface);
    ffluid_math_get_surface_spectrum(&DataCurr, &DataSpectrum);
    ffluid_math_get_volume(&DataCurr, &DataCurr.Volume);
    ffluid_math_get_hamiltonian(&DataCurr, &DataCurr.Hamiltonian, &DataCurr.SurfaceEnergy);
    ffluid_math_get_angular(&DataSurface, &DataCurr.AngularM);
    /* output with angular momentum now */
    printf("Writing %s at time %.8Qe\t", surf_name, DataCurr.time);
    printf("Vol = %.15Le\t", DataCurr.Volume);
    printf("K+P = %.15Le\t", DataCurr.Hamiltonian+DataCurr.SurfaceEnergy);
    printf("J = %.15Le\t",   DataCurr.AngularM);
    printf("z0 = %.15Le\n",  DataCurr.z0);
    /* */
    ffluid_write_surface(&DataSurface, surf_name);
    ffluid_write_spectrum(&DataSpectrum, spec_name);
    EvolveConfig.dmp_cnt++; 
  }
}

void ffluid_last_detailed_snapshot_RV() {
    //ffluid_math_get_surface_variables(&DataCurr, &DataSurface);
    sprintf(surf_name, "disc/surf_%04lu.txt", EvolveConfig.dmp_cnt);
    sprintf(full_name, "disc/full_%04lu.txt", EvolveConfig.dmp_cnt);
    sprintf(spec_name, "disc/spec_%04lu.txt", EvolveConfig.dmp_cnt);
    ffluid_math_get_surface_variables_RV(&DataCurr, &DataSurface);
    ffluid_math_get_surface_spectrum(&DataCurr, &DataSpectrum);
    ffluid_math_get_volume_RV(&DataCurr, &DataCurr.Volume);
    ffluid_math_get_hamiltonian_RV(&DataCurr, &DataCurr.Hamiltonian, &DataCurr.SurfaceEnergy);
    ffluid_math_get_angular_RV(&DataSurface, &DataCurr.AngularM);
    /* output with angular momentum */
    ffluid_append_to_log(&DataSurface, "time.log");
    printf("Writing %s at time %.5Qe\t", surf_name, DataCurr.time);
    printf("Vol = %.15Le\t", DataCurr.Volume);
    printf("K+P = %.15Le\t", DataCurr.Hamiltonian+DataCurr.SurfaceEnergy);
    printf("J = %.10Le\t", DataCurr.AngularM);
    printf("x0 = %15.8Le\t", creall(DataCurr.z0));
    printf("y0 = %15.8Le\n", cimagl(DataCurr.z0));
    ffluid_write_surface(&DataSurface, surf_name);
    ffluid_write_full_data(&DataSurface, &DataCurr, full_name);
    ffluid_write_spectrum(&DataSpectrum, spec_name);
    EvolveConfig.dmp_cnt++; 
}

void ffluid_detailed_snapshot_RV() {
  if ((EvolveConfig.cur_step % 128) == 0) {
    ffluid_math_get_surface_variables_RV(&DataCurr, &DataSurface);
    ffluid_math_get_volume_RV(&DataCurr, &DataCurr.Volume);
    ffluid_math_get_hamiltonian_RV(&DataCurr, &DataCurr.Hamiltonian, &DataCurr.SurfaceEnergy);
    ffluid_math_get_angular_RV(&DataSurface, &DataCurr.AngularM);
    ffluid_append_to_log(&DataSurface, "time.log");
  }
  if ((EvolveConfig.cur_step % 512) == 0) {
    sprintf(surf_name, "disc/surf_%04lu.txt", EvolveConfig.dmp_cnt);
    sprintf(full_name, "disc/full_%04lu.txt", EvolveConfig.dmp_cnt);
    sprintf(spec_name, "disc/spec_%04lu.txt", EvolveConfig.dmp_cnt);
    ffluid_math_get_surface_variables_RV(&DataCurr, &DataSurface);
    ffluid_math_get_surface_spectrum(&DataCurr, &DataSpectrum);
    ffluid_math_get_surface_spectrum(&DataSurface, &DataSpectrumZPh);
    ffluid_math_get_volume_RV(&DataCurr, &DataCurr.Volume);
    ffluid_math_get_hamiltonian_RV(&DataCurr, &DataCurr.Hamiltonian, &DataCurr.SurfaceEnergy);
    ffluid_math_get_angular_RV(&DataSurface, &DataCurr.AngularM);
    /* output with angular momentum now */
    printf("Writing %s at time %.5Qe\t", surf_name, DataCurr.time);
    printf("Vol = %.15Le\t", DataCurr.Volume);
    printf("K+P = %.15Le\t", DataCurr.Hamiltonian+DataCurr.SurfaceEnergy);
    printf("J = %.10Le\t",   DataCurr.AngularM);
    printf("x0 = %15.8Le\t", creall(DataCurr.z0));
    printf("y0 = %15.8Le\n", cimagl(DataCurr.z0));
    /* */
    ffluid_write_surface(&DataSurface, surf_name);
    ffluid_write_full_data(&DataSurface, &DataCurr, full_name);
    ffluid_write_full_spectrum(&DataSpectrum, &DataSpectrumZPh, spec_name);
    EvolveConfig.dmp_cnt++; 
  }
}

void ffluid_evolve() {
  ffluid_start_log("time.log");
  ffluid_detailed_snapshot_RV();  					// changed to RV 
  //ffluid_detailed_snapshot();  					
  while (EvolveConfig.cur_step < EvolveConfig.nsteps) {
    ffluid_data_copy(&DataCurr, &DataPrev);
    //ffluid_runge_kutta_4(&DataCurr);
    ffluid_runge_kutta_4_RV(&DataCurr);
    //ffluid_mapping_test_resolved(fname);
    ffluid_detailed_snapshot_RV();					// changed to RV
    //ffluid_detailed_snapshot();					
    //if (EvolveConfig.cur_step == 32768) break;
  }
  ffluid_last_detailed_snapshot_RV();
  //ffluid_last_detailed_snapshot();
  
  //ffluid_math_get_volume(&DataCurr, &DataCurr.Volume);
}
