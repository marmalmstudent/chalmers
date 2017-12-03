#include <check.h>
#include <stdlib.h>
#include "../../src/util/array.h"

static double *dcreate(double d)
{
  double *dp = malloc(sizeof(double));
  *dp = d;
  return dp;
}
static void dfree(void *d)
{
  if (d != NULL)
    free(d);
}

START_TEST(test_create)
{
  Array *arr = array_create((void (*)(void *))dfree);
  ck_assert_ptr_ne(arr, NULL);
  double *d, *dget;
  
  // add & remove entry
  d = dcreate(5);
  array_add(arr, d);
  dget = array_get(arr, 0);
  ck_assert_ptr_eq(d, dget);
  array_remove_entry(arr, d);
  dget = array_get(arr, 0);
  ck_assert_ptr_eq(NULL, dget);
  
  // add & remove idx
  d = dcreate(5);
  array_add(arr, d);
  dget = array_get(arr, 0);
  ck_assert_ptr_eq(d, dget);
  array_remove_index(arr, 0);
  dget = array_get(arr, 0);
  ck_assert_ptr_eq(NULL, dget);
  
  array_free(arr);
}
END_TEST

START_TEST(test_iterator)
{
  Array *arr = array_create(NULL);
  ck_assert_ptr_ne(arr, NULL);
  double *dget;
  
  // add & remove entry
  double d[100];
  for (unsigned i = 0; i < 100; ++i) {
    d[i] = i << 1;
    array_add(arr, &d[i]);
  }
  Iterator *itr = array_get_iterator(arr);
  for (unsigned i = 0; itr->has_next(itr); ++i) {
    dget = itr->next(itr);
    ck_assert_ptr_eq(dget, &d[i]);
  }
  itr->destroy(itr);

  array_free(arr);
}
END_TEST

extern TCase *array_case(void)
{
  TCase *tc = tcase_create("Array");
  
  tcase_add_test(tc, test_create);
  tcase_add_test(tc, test_iterator);
  
  return tc;
}
