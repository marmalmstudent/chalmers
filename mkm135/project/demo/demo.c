#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>

uint64_t systime_ms(void);
double drand(void);
void fillrand(double *arr, size_t len);

/* system time since jan 1 1970 with ms precision */
uint64_t systime_ms(void)
{
    struct timespec tstr;
    clock_gettime(CLOCK_REALTIME, &tstr);
    return (((uint64_t) tstr.tv_sec) * 1000)
         + (((uint64_t) tstr.tv_nsec) / 1000000);
}

double drand(void)
{
  return ((double) rand()) / (((double) RAND_MAX) + 1);
}

void fillrand(double *arr, size_t len)
{
  while (len--) {
    *arr++ = drand();
  }
}

void matsqrt_fast(double *arr, size_t rows, size_t cols)
{
  size_t i, j;
  for (i = 0; i < rows; ++i)
    for (j = 0; j < cols; ++j)
      *(arr + i*cols + j) = sqrt(*(arr + i*cols + j));
}

void matsqrt_slow(double *arr, size_t rows, size_t cols)
{
  size_t i, j;
  for (j = 0; j < cols; ++j)
    for (i = 0; i < rows; ++i)
      *(arr + i*cols + j) = sqrt(*(arr + i*cols + j));
}

uint64_t test(void (*f)(double *, size_t, size_t), size_t rows, size_t cols)
{
  double *arr;
  uint64_t t1, t2;
  arr = (double *) calloc(rows * cols, sizeof(double));
  fprintf(stdout, "Rows:        %10lu\n"
	          "Cols:        %10lu\n"
                  "Matrix size: %10lu MB\n",
	  rows, cols, rows*cols*sizeof(double) / (1024*1024));
  fillrand(arr, rows * cols);
  t1 = systime_ms();
  (*f)(arr, rows, cols);
  t2 = systime_ms();
  free(arr);
  return t2 - t1;
}

int main(int argc, char **argv)
{
  size_t row, col;
  
  srand(time(NULL));
  row = (size_t) atol(*(argv + 1));
  col = (size_t) atol(*(argv + 2));
  if (!strcmp(*(argv + 3), "slow")) {
    printf("Running slow code (more cache misses)\n");
    printf("Slow traverse: %" PRIu64 " ms\n", test(matsqrt_slow, row, col));
  }
  else if (!strcmp(*(argv + 3), "fast")) {
    printf("Running fast code (more cache hits)\n");
    printf("Fast traverse: %" PRIu64 " ms\n", test(matsqrt_fast, row, col));
  }
}
