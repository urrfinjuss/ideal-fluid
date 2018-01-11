
/* declared in math/equations.c */
extern void ffluid_equations_module();
extern void ffluid_call_rhs(data_ptr in, data_ptr out);

/* declared in math/halfplane_surface.c */
/* declared in math/disc_surface.c */
extern void ffluid_math_init_surface();
extern void ffluid_math_get_surface_variables(data_ptr in, data_ptr out);
extern void ffluid_math_set_zero_mode(data_ptr in, long_complex_t S0, long_complex_t *out);

