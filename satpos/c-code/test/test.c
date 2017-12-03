#include <check.h>
#include <stdlib.h>
#include "test.h"


static Suite *satpos_suite(void)
{
  Suite *s = suite_create("SatPos");

  TCase *tc_file = file_case();
  TCase *tc_array = array_case();
  TCase *tc_str = str_case();
  TCase *tc_csv_reader = csv_reader_case();
  TCase *tc_coordinate = coordinate_case();
  TCase *tc_coordinate_factory = coordinate_factory_case();
  suite_add_tcase(s, tc_file);
  suite_add_tcase(s, tc_array);
  suite_add_tcase(s, tc_str);
  suite_add_tcase(s, tc_csv_reader);
  suite_add_tcase(s, tc_coordinate);
  suite_add_tcase(s, tc_coordinate_factory);
  
  return s;
}

int main(void)
{
  Suite *suite = satpos_suite();
  SRunner *runner = srunner_create(suite);
  srunner_run_all(runner, CK_VERBOSE);
  int number_failed = srunner_ntests_failed(runner);
  srunner_free(runner);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
