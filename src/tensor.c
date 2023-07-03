#include <check.h>
#include <logging.h>
#include <tensor.h>
#include <mpi_plan.h>

MPI_Win tens_wins[MAX_WINS] = {MPI_WIN_NULL};
void *tens_adds[MAX_WINS] = {NULL};

int findAvailWin()
{
  for (int loc = 0;; loc++)
  {
    check(loc < MAX_WINS);
    if (tens_wins[loc] == MPI_WIN_NULL)
      return loc;
  }
}

int findWinByAdd(void *p)
{
  for (int loc = 0;; loc++)
  {
    check(loc < MAX_WINS);
    if (tens_adds[loc] == p)
      return loc;
  }
}

/** Returns a pointer to a 3D tensor of type dtype_sz with access range:
 * p3[st[0] - en[0]] [st[1] - en[1]] [st[2] - en[2]]
 */
void ***tensor3(int st[3], int en[3], size_t dtype_sz)
{
  int n0 = en[0] - st[0] + 1;
  int n1 = en[1] - st[1] + 1;
  int n2 = en[2] - st[2] + 1;
  void *p = fftw_malloc(n0 * n1 * n2 * dtype_sz);

  LOG_DEBUG("p = %p\n", p);

  void ***p3 = tensor3View(p, st, en, dtype_sz);
  return p3;
}

void ***tensor3Shared(int st[3], int en[3], size_t dtype_sz)
{
  int n0 = en[0] - st[0] + 1;
  int n1 = en[1] - st[1] + 1;
  int n2 = en[2] - st[2] + 1;

  int loc = findAvailWin();

  void *p;
  MPI_Win_allocate_shared(sm_rank == 0 ? n0 * n1 * n2 * dtype_sz : 0,
                          dtype_sz, MPI_INFO_NULL, sm_comm, &p, &tens_wins[loc]);

  int tmp, tmp2;
  MPI_Win_shared_query(tens_wins[loc], 0, &tmp, &tmp2, &p);
  tens_adds[loc] = p;
  LOG_DEBUG("p = %p\n", p);

  void ***p3 = tensor3View(p, st, en, dtype_sz);
  return p3;
}

/** Returns a 3D view into a block of memory pointed by t
 * p3[st[0] - en[0]] [st[1] - en[1]] [st[2] - en[2]]
 */
void ***tensor3View(void *p, int st[3], int en[3], size_t dtype_sz)
{
  int n0 = en[0] - st[0] + 1;
  int n1 = en[1] - st[1] + 1;
  int n2 = en[2] - st[2] + 1;

  /* Offset data */
  p -= st[2] * dtype_sz;

  /* dim0 x dim1 array for pointers to data */
  void **p2 = malloc(n0 * n1 * sizeof(void *));
  check(p2 != NULL);
  LOG_DEBUG("p2 = %p\n", p2);

  /* dim0 array for pointers to p2 */
  void ***p3 = malloc(n0 * sizeof(void **));
  check(p3 != NULL);
  LOG_DEBUG("p3 = %p\n", p3);

  /* match dim0 x dim1 ptrs to data */
  *p2 = p;
  for (int i = 1; i < n0 * n1; i++)
    p2[i] = p2[i - 1] + n2 * dtype_sz;
  p2 -= st[1];

  /* match dim0 pts to p2 */
  *p3 = p2;
  for (int i = 1; i < n0; i++)
    p3[i] = p3[i - 1] + n1;
  p3 -= st[0];
  // LOG_DEBUG("p3[st0] = %p\n", p3[st[0]]);

  return p3;
}

void fenceTensor3(void *p)
{
  int loc = findWinByAdd(p);
  MPI_Win_fence(0, tens_wins[loc]);
}

void freeTensor3(void ***p3, int st[3], size_t dtype_sz)
{
  p3 += st[0];
  void **p2 = *p3 + st[1];
  void *p = *p2 + st[2] * dtype_sz;
  // LOG_DEBUG("free p3 -> %p\n", p3);
  // LOG_DEBUG("free p2 -> %p\n", p2);
  // LOG_DEBUG("free p  -> %p\n", p);

  free(p3);
  free(p2);
  fftw_free(p);
}

void freeTensor3View(void ***p3, int st[3], size_t dtype_sz)
{
  p3 += st[0];
  void **p2 = *p3 + st[1];
  void *p = *p2 + st[2] * dtype_sz;
  // LOG_DEBUG("free p3 -> %p\n", p3);
  // LOG_DEBUG("free p2 -> %p\n", p2);
  // LOG_DEBUG("free p  -> %p\n", p);

  free(p3);
  free(p2);
}

void freeTensor3Shared(void ***p3, int st[3], size_t dtype_sz)
{
  p3 += st[0];
  void **p2 = *p3 + st[1];
  void *p = *p2 + st[2] * dtype_sz;

  int loc = findWinByAdd(p);

  free(p3);
  free(p2);
  MPI_Win_free(&tens_wins[loc]);
}
