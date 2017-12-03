#include "array.h"
#include <stdlib.h>
#include <stdio.h>

struct array_t {
  size_t capacity;
  size_t length;
  void **entry;
  void (*delete)(void *);
};

struct iterator_data_t {
  size_t offset;
};

enum { INITIAL_CAPACITY=8 };

static int array_iterator_hasnext(Iterator *itr);
static void *array_iterator_next(Iterator *itr);
static void array_iterator_remove(Iterator *itr);

static void increase_capacity(Array *arr);

extern Array *array_create(void (*delete)(void *))
{
  Array *root = malloc(sizeof(Array));
  root->capacity = INITIAL_CAPACITY;
  root->entry = calloc(root->capacity, sizeof(void *));
  root->length = 0;
  root->delete = delete;
  return root;
}

extern void array_free(Array *arr)
{
  if (arr->entry != NULL) {
    for (size_t i = 0; i < arr->capacity; ++i) {
      if (arr->entry[i] != NULL) {
	if (arr->delete != NULL)
	  arr->delete(arr->entry[i]);
	arr->entry[i] = NULL;
      }
    }
    free(arr->entry);
    arr->entry = NULL;
  }
  free(arr);
}

extern void array_add(Array *arr, void *entry)
{
  if (arr->length >= arr->capacity)
    increase_capacity(arr);
  arr->entry[arr->length++] = entry;
}

extern void array_remove_entry(Array *arr, const void *entry)
{
  for (size_t i = 0; i < arr->length; ++i)
    if (arr->entry[i] == entry) {
      array_remove_index(arr, i);
      break;
    }
}

extern void array_remove_index(Array *arr, size_t index)
{
  if (index < arr->length) {
      arr->delete(arr->entry[index]);
      size_t i = index;
      for(; i < arr->length-1; ++i)
	arr->entry[i] = arr->entry[i+1];
      arr->entry[i] = NULL;
      arr->length--;
  }
}
extern void *array_get(Array *arr, size_t index)
{
  return index < arr->length ? arr->entry[index] : NULL;
}

extern Iterator *array_get_iterator(Array *host)
{
  Iterator *itr = malloc(sizeof(Iterator));
  itr->data = malloc(sizeof(struct iterator_data_t));
  itr->data->offset = 0;
  itr->destroy = array_iterator_free;
  
  itr->host = host;
  itr->has_next = array_iterator_hasnext;
  itr->next = array_iterator_next;
  itr->remove = array_iterator_remove;
  return itr;
}

extern void array_iterator_free(Iterator *itr)
{
  if (itr->data != NULL)
    free(itr->data);
  free(itr);
}

static int array_iterator_hasnext(Iterator *itr)
{
  return itr->data->offset < ((Array *) itr->host)->length;
}

static void *array_iterator_next(Iterator *itr)
{
  return itr->has_next(itr) ? ((Array *) itr->host)->entry[itr->data->offset++] : NULL;
}

static void array_iterator_remove(Iterator *itr)
{
  array_remove_index(itr->host, itr->data->offset);
}

static void increase_capacity(Array *arr)
{
  size_t capacity = arr->capacity;
  arr->capacity = arr->capacity + (arr->capacity >> 1);
  arr->entry = realloc(arr->entry, arr->capacity*sizeof(void *));
  for (size_t i = capacity; i < arr->capacity; ++i)
    arr->entry[i] = NULL;
}
