/* declared in io/io.c */
extern void ffluid_io_module();

/* declared in io/input_data.c */
extern void ffluid_read_initial_data(data_ptr in);
extern void ffluid_set_initial_data(data_ptr in);
extern void ffluid_read_mapping_parameters(data_ptr in, FILE *fh);
extern void ffluid_read_data_from_file(data_ptr in, FILE *fh);

/* declared in io/input_parameters.c */
extern void ffluid_read_cl_arguments(int narg, char **argv);  
extern void ffluid_read_input_file(char *filename);
extern void ffluid_scan_input_file(FILE *fh);

/* declared in io/output.c */
extern void ffluid_write_surface(data_ptr in, char *fname);
extern void ffluid_write_spectrum(data_ptr in, char *fname);
extern void ffluid_start_log(char *fname);
extern void ffluid_append_to_log(data_ptr in, char *fname);
extern void ffluid_write_array(long_complex_t *in, unsigned long N, char *fname);
