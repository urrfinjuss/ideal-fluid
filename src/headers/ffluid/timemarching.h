/* declared in timemarching/stepping.c */
extern void ffluid_timemarching_module();
extern void ffluid_setup_stepping();
extern void ffluid_evolve();

/* declared in timemarching/runge_kutta_4.c */
extern void ffluid_init_runge_kutta_4();
extern void ffluid_runge_kutta_4(data_ptr in);
extern void ffluid_runge_kutta_4_RV(data_ptr in);

/* declared in timemarching/runge_kutta_6.c */
extern void ffluid_init_runge_kutta_6();
extern void ffluid_runge_kutta_6(data_ptr in);

/* Shared variables for time-marching subsystem */
extern data_ptr tmp;
extern data_ptr *rhs;

