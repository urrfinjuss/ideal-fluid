#include "config.h"
#ifdef USE_QUAD
#include "ffluidq.h"
//#include "ffluidq/ffluidq_mapping.h"
#else
#include "ffluid.h"
//#include "ffluid/ffluid_mapping.h"
#endif

sim_data Data;
evolve_params Evolve_Config;

int main (int argc, char **argv) {
  data_ptr Data_Ptr = &Data; 
  ffluid_read_cl_arguments(argc, argv);
  
  ffluid_disk();
  ffluid_halfplane();
  
  ffluid_init_data(Data_Ptr);
  return 0;
}
