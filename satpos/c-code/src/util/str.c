#include "str.h"
#include <string.h>
#include <stdlib.h>


static void dfree(char *c)
{
  if (c != NULL)
    free(c);
}

static void add_string(Array *arr, const char *text, size_t len)
{
  char *str = calloc(len + 1, sizeof(char));
  strncpy(str, text, len);
  str[len] = '\0';
  array_add(arr, str);
}

extern Array *str_split(const char *text, const char *delim)
{
  Array *arr = array_create((void (*)(void *))dfree);
  size_t len = strlen(text);
  size_t delim_len = strlen(delim);
  if (delim_len > len)
    return arr;
  size_t i = 0, j = 0;
  while (i <= len - delim_len) {
    if (!strncmp(text + i, delim, delim_len)) {
      add_string(arr, text+j, i-j);
      i += delim_len;
      j = i;
    } else {
      ++i;
    }
  }
  add_string(arr, text+j, i-j);
  return arr;
}
