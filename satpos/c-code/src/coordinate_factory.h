#ifndef COORDINATE_FACTORY_H
#define COORDINATE_FACTORY_H

typedef struct __coordinate_factory_t CoordinateFactory;

extern CoordinateFactory *coordinate_factory_create(void);
extern void coordinate_factory_free(CoordinateFactory *cf);

#endif //COORDINATE_FACTORY_H
