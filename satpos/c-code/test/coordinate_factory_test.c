#include <check.h>
#include <stdlib.h>
#include <stdio.h>
#include "../src/coordinate_factory.h"
#include "../src/util/str.h"


START_TEST(test_create)
{
  CoordinateFactory *fc = coordinate_factory_create();
  ck_assert_ptr_ne(fc, NULL);
  coordinate_factory_free(fc);
}
END_TEST

extern TCase *coordinate_factory_case(void)
{
  TCase *tc = tcase_create("CoordinateFactory");
  
  tcase_add_test(tc, test_create);
  
  return tc;
}
