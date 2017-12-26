typedef struct stepping_parameters {
  unsigned long nsteps;
  double 	cfl, dt;
  double	final_time;
} evolve_params, *evolve_params_ptr;

/* Shared with time marching subsystem */
extern evolve_params Evolve_Config;

extern void init_stepping(double final_time, double max_cfl);

