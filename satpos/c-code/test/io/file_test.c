#include <check.h>
#include <stdlib.h>
#include "../../src/io/file.h"

START_TEST(test_create)
{
  File *f;
  f = file_create("testfile.txt");
  ck_assert_str_eq(f->filepath, "testfile.txt");
  ck_assert_ptr_ne(f->buffer, NULL);
  ck_assert_ptr_eq(f->contents, NULL);
  ck_assert_ptr_eq(f->fp, NULL);
  
  file_free(f);
}
END_TEST

START_TEST(test_read)
{
  File *f;
  f = file_create("testfile.txt");
  read(f);
  ck_assert_str_eq(f->contents, "file contents\n");
  
  file_free(f);
}
END_TEST

extern TCase *file_case(void)
{
  TCase *tc = tcase_create("File");
  
  tcase_add_test(tc, test_create);
  tcase_add_test(tc, test_read);
  
  return tc;
}
