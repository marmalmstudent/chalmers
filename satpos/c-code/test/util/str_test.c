#include <check.h>
#include <stdlib.h>
#include <stdio.h>
#include "../../src/util/str.h"


START_TEST(test_split_single_delim)
{
  char *c = ",1.2345,2.3456,3.4567,SASK,";
  char *carr[] = {"","1.2345","2.3456","3.4567","SASK",""};
  Array *arr = str_split(c, ",");

  char *str;
  Iterator *itr = array_get_iterator(arr);

  for (unsigned i = 0; itr->has_next(itr); ++i) {
    str = itr->next(itr);
    ck_assert_ptr_ne(str, NULL);
    ck_assert_str_eq(str, carr[i]);
  }
  itr->destroy(itr);

  array_free(arr);
}
END_TEST

START_TEST(test_split_multiple_delim)
{
  char *c = "fishThiSfishISfishAfishFishfish";
  char *carr[] = {"","ThiS","IS","A","Fish",""};
  Array *arr = str_split(c, "fish");
  
  char *str;
  Iterator *itr = array_get_iterator(arr);
  for (unsigned i = 0; itr->has_next(itr); ++i) {
    str = itr->next(itr);
    ck_assert_ptr_ne(str, NULL);
    ck_assert_str_eq(str, carr[i]);
  }
  itr->destroy(itr);
  
  array_free(arr);
}
END_TEST

extern TCase *str_case(void)
{
  TCase *tc = tcase_create("Str");
  
  tcase_add_test(tc, test_split_single_delim);
  tcase_add_test(tc, test_split_multiple_delim);
  
  return tc;
}
