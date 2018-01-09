/* declared in memory/memory.c */
extern void ffluid_memory_module();
extern void ffluid_init_data(data_ptr in);
extern void ffluid_alloc_aux_array(aux_data_ptr in, unsigned long NArrays, unsigned long NElements);
extern void ffluid_alloc_fft_plans(aux_data_ptr in, fft_list_ptr out);
extern void ffluid_dealloc_aux_array(aux_data_ptr in);
