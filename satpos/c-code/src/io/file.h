#ifndef __FILE_H
#define __FILE_H

#include <stdio.h>

typedef struct __file_t File;

struct __file_t {
  FILE *fp;
  char *filepath;
  char *contents;
  size_t file_size;
  struct __buffer_t *buffer;
};

extern File *file_create(const char *filepath);
extern void file_free(File *f);

/**
 * places the contents of a file in f.
 * opentype: "r", "r+", "w", "w+", "a", "a+"
 */
extern int read(File *f);

#endif // __FILE_H
