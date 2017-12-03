#include <check.h>
#include <stdlib.h>
#include <stdio.h>
#include "../../src/util/str.h"
#include "../../src/io/csv_reader.h"


START_TEST(test_read_csv)
{
  Array *rows = read_csv("../../res/csv/BAKE.lon.csv");
  /*
  Iterator *itr1 = array_get_iterator(rows);
  while (itr1->has_next(itr1)) {
    Iterator *itr2 = array_get_iterator(itr1->next(itr1));
    printf("[ ");
    while (itr2->has_next(itr2)) {
      printf("%6.2f ", *((double *)itr2->next(itr2)));
    }
    printf("]\n");
    itr2->destroy(itr2);
  }
  itr1->destroy(itr1);
  array_free(rows);
  */
}
END_TEST

extern TCase *csv_reader_case(void)
{
  TCase *tc = tcase_create("CSVReader");
  
  tcase_add_test(tc, test_read_csv);
  
  return tc;
}
