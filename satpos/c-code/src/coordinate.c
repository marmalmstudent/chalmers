#include "coordinate.h"
#include <stdlib.h>

static Point *point_create(void);
static void point_free(Point *p);

extern Coordinate *coordinate_create(void)
{
  Coordinate *c = (Coordinate *) malloc(sizeof(Coordinate));
  c->time = NULL;
  c->lon = NULL;
  c->lat = NULL;
  c->rad = NULL;
  return c;
}

extern void coordinate_free(Coordinate *c)
{
  if (c->time != NULL) { free(c->time); c->time = NULL; }
  if (c->lon != NULL) { point_free(c->lon); c->lon = NULL; }
  if (c->lat != NULL) { point_free(c->lat); c->lat = NULL; }
  if (c->rad != NULL) { point_free(c->rad); c->rad = NULL; }
  free(c);
}

static Point *point_create(void)
{
  Point *p = (Point *) malloc(sizeof(Point));
  return p;
}

static void point_free(Point *p)
{
  free(p);
}
