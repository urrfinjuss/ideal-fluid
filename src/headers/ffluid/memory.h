/* declared in memory/memory.c */
extern void ffluid_memory_module();
extern void ffluid_init_data();
extern void ffluid_init_grid(data_ptr in, grid_ptr out);
extern void ffluid_reinit_grid(data_ptr in, grid_ptr out);
/* undeclared */
extern void allocate_data(data_ptr in);
extern void deallocate_data(data_ptr in);
