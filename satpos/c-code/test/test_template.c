#include <check.h>
#include <stdlib.h>

START_TEST(test_addition)
{
  ck_assert_int_eq(1+1, 2);
}
END_TEST

extern TCase *template_case(void)
{
  TCase *tc = tcase_create("Template");
  
  tcase_add_test(tc, test_addition);
  
  return tc;
}
