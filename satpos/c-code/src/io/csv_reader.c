#include "csv_reader.h"
#include <string.h>
#include <stdlib.h>
#include "../util/str.h"
#include "file.h"


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

extern Array *read_csv(const char *fpath)
{
  File *f = file_create(fpath);
  read(f);
  
  Array *lines = str_split(f->contents, "\n");
  Array *rows = array_create((void (*)(void *))array_free);
  Iterator *itr1 = array_get_iterator(lines);
  while (itr1->has_next(itr1)) {
    const char *linestr = itr1->next(itr1);
    if (strlen(linestr) > 0) {
      Array *line = str_split(linestr, ",");
      Array *cols = array_create(dfree);
      Iterator *itr2 = array_get_iterator(line);
      for (unsigned i = 0; i < 3 && itr2->has_next(itr2); ++i)
	array_add(cols, dcreate(strtod(itr2->next(itr2), NULL)));
      itr2->destroy(itr2);
      array_free(line);
      array_add(rows, cols);
    }
  }
  itr1->destroy(itr1);
  printf("Freeing lines: %p\n", lines);
  array_free(lines);
  printf("Freed lines: %p\n", lines);

  file_free(f);
  return rows;
}
