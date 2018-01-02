/* declared in io/io.c */
extern void ffluid_io_module();

/* declared in io/input.c */
extern void ffluid_scan_input_file(FILE *fh);
extern void ffluid_read_input_file(char *filename);
extern void ffluid_read_cl_arguments(int narg, char **argv);  
extern void ffluid_read_initial_data();
extern void ffluid_set_initial_data();
extern void ffluid_read_mapping_parameters(FILE *fh);
extern void ffluid_read_data_from_file(FILE *fh);

/* declared in io/output.c */
extern void ffluid_write_surface(data_ptr in, char *fname);
extern void ffluid_append_to_log();
