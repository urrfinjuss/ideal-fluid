#include <config.h>
#ifdef USE_QUAD
#include "fffluidq.h"
#else
#include "fffluid.h"
#endif

int main (void)
{
  printf("This is %s\n", PACKAGE_STRING);
  fffluid_disk();
  fffluid_halfplane();
  return 0;
}
