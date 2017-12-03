#include <check.h>
#include <stdlib.h>
#include "../src/coordinate.h"

START_TEST(test_create)
{
  Coordinate *c = coordinate_create();
  ck_assert_ptr_eq(c->time, NULL);
  ck_assert_ptr_eq(c->lon, NULL);
  ck_assert_ptr_eq(c->lat, NULL);
  ck_assert_ptr_eq(c->rad, NULL);

  coordinate_free(c);
}
END_TEST

extern TCase *coordinate_case(void)
{
  TCase *tc = tcase_create("Coordinate");
  
  tcase_add_test(tc, test_create);
  
  return tc;
}
