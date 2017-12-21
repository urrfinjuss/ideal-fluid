#include <config.h>
#include "src/header.h"
#include "src/mapping/mapping.h"

int main (void)
{
  printf("This is %s\n", PACKAGE_STRING);
  dummy_disk();
  dummy_halfplane();
  return 0;
}
