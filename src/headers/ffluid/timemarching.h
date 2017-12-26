extern void init_runge_kutta_4();
extern void init_runge_kutta_6();
extern void runge_kutta_4(data_ptr in, const double step_size);
extern void runge_kutta_6(data_ptr in, const double step_size);
extern void evolve(data_ptr in, const double final_time);

/* Shared variables for time-marching subsystem */
extern data_ptr tmp;
extern data_ptr *rhs;

