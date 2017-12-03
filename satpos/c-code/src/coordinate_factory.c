#include "coordinate_factory.h"
#include <stdlib.h>

struct __coordinate_factory_t {
  double *lat;
  double *lon;
  double *rad;
};

extern CoordinateFactory *coordinate_factory_create(void)
{
  CoordinateFactory *cf = (CoordinateFactory *) malloc(sizeof(CoordinateFactory));
  return cf;
}

extern void coordinate_factory_free(CoordinateFactory *cf)
{
  free(cf);
}
