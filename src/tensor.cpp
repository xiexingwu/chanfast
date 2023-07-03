#include <mpi.h>
#include <fftw3.h> // Maybe the wrapper fftw++.h?

#include <spdlog/spdlog.h>

#include <chanfast/check.h>
#include <chanfast/tensor.h>
#include <chanfast/mpi_plan.h>

#define MAX_WINS 64

MPI_Win tens_wins[MAX_WINS] = {MPI_WIN_NULL};
int findAvailWin();

// Tensor3
template <typename T>
Tensor3<T>::Tensor3() {} // Dummy

template <typename T>
Tensor3<T>::Tensor3(
    int st[3], int en[3], bool shared)
{
  this->shared = shared;

  int n0 = en[0] - st[0] + 1;
  int n1 = en[1] - st[1] + 1;
  int n2 = en[2] - st[2] + 1;
  int ncells = n0 * n1 * n2;

  T *p;

  if (shared)
  {
    id = findAvailWin();

    MPI_Win_allocate_shared(Mpi::sm_rank == 0 ? ncells * sizeof(T) : 0,
                            sizeof(T), MPI_INFO_NULL, Mpi::sm_cm, &p, &tens_wins[id]);

    int tmp, tmp2;
    MPI_Win_shared_query(tens_wins[id], 0, &tmp, &tmp2, &p);
  }
  else
  {
    p = fftw_malloc(ncells * sizeof(T));
  }

  view = View3(p, st, en);
}

template <typename T>
Tensor3<T>::~Tensor3()
{
  if (shared)
  {
    MPI_Win_free(&tens_wins[id]);
  }
  {
    fftw_free(view.p);
  }

  delete view;
}

template <typename T>
template <typename U>
View3<T>::View3(
    const U *p,
    int st[3], int en[3])
{
  int n0 = en[0] - st[0] + 1;
  int n1 = en[1] - st[1] + 1;
  int n2 = en[2] - st[2] + 1;
  // TODO - This doens't feel right.
  // pointers st en get moved, but still got allocated space for 3 ints on allocating View?
  this->st = st;
  this->en = en;
  this->ncells = n0 * n1 * n2;

  this->p = (T *)p;

  // dim0 x dim1 array for pointers to data
  p2 = new T *[n0 * n1];
  *p2 = (T *)p - st[2]; // offset z
  for (int i = 1; i < n0 * n1; i++)
    p2[i] = p2[i - 1] + n2;

  // dim0 array for pointers to p2
  p3 = new T **[n0];
  *p3 = p2 - st[1]; // offset y
  for (int i = 1; i < n0; i++)
    p3[i] = p3[i - 1] + n1;

  // Only p3 is manually offset from allocation since it is used in 3d global indexing.
  p3 -= st[0]; // offset x
}
template <typename T>
View3<T>::~View3()
{
  p3 += st[0];

  delete[] p3;
  delete[] p2;
}
// Utility
int findAvailWin()
{
  for (int loc = 0;; loc++)
  {
    check(loc < MAX_WINS);
    if (tens_wins[loc] == MPI_WIN_NULL)
      return loc;
  }
}
