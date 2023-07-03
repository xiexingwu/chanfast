#include <check.h>
#include <logging.h>

#include <stdarg.h>

static MPI_File fp = MPI_FILE_NULL;
static char buf[512];

void initLogging()
{
  check(RANK != -1);
  char fname[256];
  sprintf(fname, "logs/%d.log", RANK);
  MPI_File_open(MPI_COMM_SELF, fname,
    MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_UNIQUE_OPEN,
    MPI_INFO_NULL, &fp);
  check(fp != MPI_FILE_NULL);
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
    fflush(stderr);
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
  check(fp != MPI_FILE_NULL);
  int cnt = sprintf(buf, "[INFO] ");
  va_list args;
  va_start(args, fmt);
  cnt += vsprintf(buf+cnt, fmt, args);
  va_end(args);

  MPI_File_write(fp, buf, cnt, MPI_CHAR, MPI_STATUS_IGNORE);
}

void LOG_DEBUG(const char *fmt, ...)
{
#ifdef DEBUG
  check(fp != MPI_FILE_NULL);
  int cnt = sprintf(buf, "[DEBUG] ");
  va_list args;
  va_start(args, fmt);
  cnt += vsprintf(buf + cnt, fmt, args);
  va_end(args);
  
  MPI_File_write(fp, buf, cnt, MPI_CHAR, MPI_STATUS_IGNORE);
#endif /* DEBUG */
}

extern void freeLogging()
{
  LOG_DEBUG("Free logging...\n");
  check(fp != MPI_FILE_NULL);
  MPI_File_close(&fp);
}
