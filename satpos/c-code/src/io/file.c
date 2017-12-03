#include "file.h"
#include <stddef.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define buffsize 4096

typedef struct __buffer_t Buffer;

struct __buffer_t {
  size_t len;
  size_t offset;
  char *buffer;
};

static Buffer *buffer_create(void);
static void buffer_free(Buffer *b);
static void buffer_extend(Buffer *b);
static int open(File *f, const char *access_type);
static void close(File *f);

extern File *file_create(const char *filepath)
{
  File *f = calloc(1, sizeof(File));
  f->filepath = calloc(strlen(filepath) + 1, sizeof(char));
  strcpy(f->filepath, filepath);
  f->buffer = buffer_create();
  return f;
}

extern void file_free(File *f)
{
  if (f->buffer != NULL) {
    buffer_free(f->buffer);
    f->buffer = NULL;
  }
  f->contents = NULL;
  if (f->fp != NULL)
    close(f);
  if (f->filepath != NULL) {
    free(f->filepath);
    f->filepath = NULL;
  }
  free(f);
}

extern int read(File *f)
{
  if (open(f, "r")) {
    Buffer *b = f->buffer;
    for (int c; (c = fgetc(f->fp)) != EOF;) {
      if (b->offset >= b->len)
	buffer_extend(b);
      b->buffer[b->offset++] = (char) c;
    }
    b->buffer[b->offset] = '\0';
    f->contents = b->buffer;
    f->file_size = b->offset;
    
    close(f);
    return 1;
  }
  return 0;
}

static Buffer *buffer_create(void)
{
  Buffer *b = calloc(1, sizeof(Buffer));
  b->len = buffsize;
  buffer_extend(b);
  return b;
}

static void buffer_free(Buffer *b)
{
  if (b->buffer != NULL) {
    free(b->buffer);
    b->buffer = NULL;
  }
  free(b);
}

static void buffer_extend(Buffer *b)
{
  b->len = 3*b->len / 2;
  b->buffer = (char *) realloc(b->buffer, b->len * sizeof(char));
}

static int open(File *f, const char *access_type)
{
  if (f->fp != NULL)
    close(f);
  if ((f->fp = fopen(f->filepath, access_type)) == NULL)
      perror("Error opening file");
  return f->fp != NULL;
}

static void close(File *f)
{
  if (f->fp != NULL) {
    if (fclose(f->fp))
      perror("Error closing file");
    else
      f->fp = NULL;
  }
}
