#ifndef COORDINATE_H
#define COORDINATE_H

typedef struct __point_t Point;
typedef struct __coordinate_t Coordinate;

struct __point_t {
  double time;
  double val;
  double error;
};

struct __coordinate_t {
  double *time;
  Point *lat;
  Point *lon;
  Point *rad;
};

extern Coordinate *coordinate_create(void);
extern void coordinate_free(Coordinate *c);

#endif //COORDINATE_H
