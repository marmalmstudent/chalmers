#ifndef ITERATOR_H
#define ITERATOR_H

typedef struct iterator_t Iterator;
typedef struct iterator_data_t IteratorData;

struct iterator_t {
  IteratorData *data;
  void (*destroy)(Iterator *itr);

  void *host;
  void *(*next)(Iterator *itr);
  int (*has_next)(Iterator *itr);
  void (*remove)(Iterator *itr);
};

#endif //ITERATOR_H
