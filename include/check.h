#ifndef CHECK_H
#define CHECK_H

#include <mpi.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

extern int RANK;

#define DEBUG_PRINT0 \
  if (RANK == 0)    \
  DEBUG_PRINT

#ifdef DEBUG
#define DEBUG_PRINT(...) fprintf(stderr, __VA_ARGS__); fflush(stderr)
#else
#define DEBUG_PRINT(...) do {} while (0)
#endif


#define check(e) \
  ((void)((e) ? 0 : __check(#e, __FILE__, __LINE__)))

#define __check(e, file, line) \
  ( \
    (void)fprintf(stderr, "%s:%u: failed check `%s'\n", file, line, e), \
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE), \
    0 \
  )

#define unused(e) ((void)(e))

#define MAX(a, b) \
  ( \
    { __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
      _a > _b ? _a : _b; } \
  )

#define MIN(a, b) \
  ( \
      { __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
      _a < _b ? _a : _b; } \
  )

#endif /* CHECK_H */