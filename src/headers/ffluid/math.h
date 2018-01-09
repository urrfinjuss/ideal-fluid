
/* declared in math/equations.c */
extern void ffluid_equations_module();
extern void ffluid_call_rhs(data_ptr in, data_ptr out);

/* declared in math/halfplane_surface.c */
extern void ffluid_math_init_surface_halfplane();
extern void ffluid_math_get_surface_variables_halfplane(data_ptr in, data_ptr out);
extern void ffluid_math_set_zero_mode(data_ptr in, long_complex_t S0, long_complex_t *out);

/* decalred in math/disc_surface.c */
extern void ffluid_math_init_surface_disc();
extern void ffluid_math_get_surface_variables_disc(data_ptr in, data_ptr out);
