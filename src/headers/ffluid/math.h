
/* declared in math/equations.c */
extern void ffluid_alloc_equations();
extern void ffluid_equations_module();
extern void ffluid_call_rhs(data_ptr in, data_ptr out);
extern void ffluid_call_rhsX(data_ptr in, data_ptr out); // experimental
extern void ffluid_call_rhsXX(data_ptr in, data_ptr out); // experimental by text
extern void ffluid_call_rhsRV(data_ptr in, data_ptr out); // experimental by main_v06.pdf
extern void ffluid_filter_high(data_ptr in, const int fflag);

/* declared in math/halfplane_surface.c */
/* declared in math/disc_surface.c */
extern void ffluid_math_init_surface();
extern void ffluid_math_get_surface_variables(data_ptr in, data_ptr out);
extern void ffluid_math_set_zero_mode(data_ptr in, long_complex_t S0, long_complex_t *out);

/* declared in math/disc_surface.c */
/* deprecated */
extern void ffluid_math_get_volume(data_ptr in, long_double_t *volume);
extern void ffluid_math_get_angular(data_ptr in, long_double_t *angular_m); 
extern void ffluid_math_get_hamiltonian(data_ptr in, long_double_t *hamiltonian, long_double_t *surf_energy);
extern void ffluid_math_get_r0(data_ptr in);

/* declared in math/disc_surface.c */
extern void ffluid_math_get_volume_RV(data_ptr in, long_double_t *volume);
extern void ffluid_math_get_angular_RV(data_ptr in, long_double_t *angular_m); 
extern void ffluid_math_get_hamiltonian_RV(data_ptr in, long_double_t *hamiltonian, long_double_t *surf_energy);
extern void ffluid_math_get_surface_variables_RV(data_ptr in, data_ptr out);
extern void ffluid_math_get_surface_spectrum(data_ptr in, data_ptr out);
