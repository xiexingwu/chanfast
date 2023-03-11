#include <check.h>
#include <logging.h>

#include <stdarg.h>

static FILE *fp = NULL;
void initLogging()
{
  check(RANK != -1);
  char fname[256];
  sprintf(fname, "logs/%d.log", RANK);
  fp = fopen(fname, "w");
}

void LOG_STDOUT(int rank, const char *fmt, ...)
{
  if (rank == RANK)
  {
    printf("[%d] ", rank);
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
    fflush(stdout);
  }
}

void LOG_WARN(int rank, const char *fmt, ...)
{
  if (rank == RANK)
  {
    printf("[%d-WARN] ", rank);
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
    fflush(stdout);
  }
}

void LOG_ERR(int rank, const char *fmt, ...)
{
  if (rank == RANK)
  {
    fprintf(stderr, "[%d-ERR] ", rank);
    va_list args;
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fflush(stderr);
  }
}

void LOG_INFO(const char *fmt, ...)
{
  check(fp != NULL);
  fprintf(fp, "[INFO] ");
  va_list args;
  va_start(args, fmt);
  vfprintf(fp, fmt, args);
  va_end(args);
  fflush(fp);
}

void LOG_DEBUG(const char *fmt, ...)
{
#ifdef DEBUG
  check(fp != NULL);
  fprintf(fp, "[DEBUG] ");
  va_list args;
  va_start(args, fmt);
  vfprintf(fp, fmt, args);
  va_end(args);
  fflush(fp);
#endif /* DEBUG */
}

extern void freeLogging()
{
  check(fp != NULL);
  fclose(fp);
}
