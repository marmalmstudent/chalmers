#ifndef ARRAY_H
#define ARRAY_H

#include "iterator.h"
#include <stddef.h>

typedef struct array_t Array;

extern Array *array_create(void (*delete)(void *));
extern void array_free(Array *arr);
extern void array_add(Array *arr, void *entry);
extern void array_remove_entry(Array *arr, const void *entry);
extern void array_remove_index(Array *arr, size_t index);
extern void *array_get(Array *arr, size_t index);

extern Iterator *array_get_iterator(Array *host);
extern void array_iterator_free(Iterator *itr);

#endif //ARRAY_H
