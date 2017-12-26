/* This is the main header file for idealfluid program running in double precision */
/* for quad-precision __float128 default data-type ffluidq.h is included instead  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

#define PI acos(-1.0)

typedef struct data_array {
  unsigned long	N;
  fftw_complex	*Q, *V;
  double	q0, u0, l;
  double 	time;
} sim_data, *data_ptr;

typedef struct aux_array {
} aux_data, *aux_data_ptr;

typedef struct fft_array {
} fft_list, *fft_list_ptr;

/* Note to self: remember top to bottom this time */
#include "ffluid/memory.h"
#include "ffluid/mapping.h"
#include "ffluid/timemarching.h"
#include "ffluid/stepping.h"
#include "ffluid/array_func.h"
#include "ffluid/equations.h"
#include "ffluid/io.h"
#include "ffluid/messages.h"
/* Global Variable */ 
extern sim_data Data;
