/* This is the main header file for idealfluid program running in double precision */
/* for quad-precision __float128 default data-type ffluidq.h is included instead  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <quadmath.h>
#include <complex.h>


#ifdef USE_QUAD
  #include <fftw3q.h>
  typedef __float128	long_double_t;
  typedef fftwq_complex	long_complex_t;
  #define PI 	acosq(-1.0)
#else
  #include <math.h>
  #include <fftw3.h>
  typedef long double	long_double_t;
  typedef fftwl_complex	long_complex_t;
  #define PI 	acosl(-1.0)
#endif

typedef struct stepping_parameters {
  unsigned long nsteps, cur_step, dmp_cnt;
  __float128 	cfl, max_cfl, dt;
  __float128	final_time;
} evolve_params, *evolve_params_ptr;

typedef struct data_array {
  unsigned long		N;
  long_complex_t	*R, *V;
  long_complex_t	*u, *du;
  long_double_t		*q;
  long_double_t		q0, u0, l;
  long_double_t		Volume, Hamiltonian;
  __float128 		time;
} sim_data, *data_ptr;

#ifdef USE_QUAD
  #include "ffluidq/io.h"
#else
  #include "ffluid/io.h"
#endif

typedef struct aux_array {
  unsigned long		NElements;
  unsigned long		NArrays;
  long_complex_t	**X;
  long_complex_t	**Y;
} aux_data, *aux_data_ptr;

typedef struct fft_array {
  unsigned long		NFFTs;
  fftwl_plan		*fp;
  fftwl_plan		*bp;
} fft_list, *fft_list_ptr;

/* Note to self: remember top to bottom this time */
#include "ffluid/memory.h"
#include "ffluid/mapping.h"
#include "ffluid/math.h"
#include "ffluid/timemarching.h"
#include "ffluid/array_func.h"
#include "ffluid/messages.h"

typedef struct control_parameters_array{
  data_ptr		DataPtrCurr, DataPtrPrev;
  evolve_params_ptr 	EvolvePtr;
  char 			run_name[80];
  char			res_name[80];
  char			data_path[80];
} control_params, *control_params_ptr;

/* Global Variable */ 
extern sim_data		DataCurr, DataPrev;
extern sim_data         DataSpectrum, DataSurface;
extern control_params	Control;
extern evolve_params	EvolveConfig;

