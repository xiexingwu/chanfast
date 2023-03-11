#include <check.h>
#include <tensor.h>
#include <mpi_plan.h>

static const int st[3] = {1,3,5};
static const int en[3] = {3,6,9};
static const int st2[3] = {1, 2, 3};
static const int en2[3] = {5, 5, 5};
static const int n0 = en[0] - st[0] + 1;
static const int n1 = en[1] - st[1] + 1;
static const int n2 = en[2] - st[2] + 1;

void testTensor3()
{
  /* Only test for srank 0 on shared-memory comm */
  if (sm_rank != 0) return;

  int ***p3 = (int ***) tensor3(st, en, sizeof(int));

  MPI_Barrier( MPI_COMM_WORLD);
  int ***q3 = viewTensor3(&p3[st[0]][st[1]][st[2]], st2, en2, sizeof(int));

  /* Write 3D */
  for (int i = st[0]; i <= en[0]; i++)
    for (int j = st[1]; j <= en[1]; j++)
      for (int k = st[2]; k <= en[2]; k++)
        p3[i][j][k] = 100 * i + 10 * j + k;

  /* Read 1D */
  for (int i = 0; i < n0 * n1 * n2; i++)
    LOG_DEBUG("p[%2d] (%p) = %d\n", i, &p3[st[0]][st[1]][st[2] + i], p3[st[0]][st[1]][st[2] + i]);

  /* Write 3D (view) */
  for (int i = st2[0]; i <= en2[0]; i++)
    for (int j = st2[1]; j <= en2[1]; j++)
      for (int k = st2[2]; k <= en2[2]; k++)
        q3[i][j][k] = 100 * i + 10 * j + k;

  /* Read 1D */
  for (int i = 0; i < n0 * n1 * n2; i++)
    LOG_DEBUG("q[%2d] (%p) = %d\n", i, &p3[st[0]][st[1]][st[2] + i], p3[st[0]][st[1]][st[2] + i]);
  freeTensor3(p3, st, sizeof(int));
}


void testTensor3Shared(){
  int ***p3 = (int ***)tensor3Shared(st, en, sizeof(int));

  int ***q3 = viewTensor3(&p3[st[0]][st[1]][st[2]], st2, en2, sizeof(int));

  /* Write 3D */
  if (sm_rank == 0)
    for (int i = st[0]; i <= en[0]; i++)
      for (int j = st[1]; j <= en[1]; j++)
        for (int k = st[2]; k <= en[2]; k++)
          p3[i][j][k] = 100 * i + 10 * j + k;

  fenceTensor3(&p3[st[0]][st[1]][st[2]]);
  LOG_DEBUG("Finished write.\n");
  MPI_Barrier(MPI_COMM_WORLD);
  /* Read 1D */
  for (int i = 0; i < n0 * n1 * n2; i++)
    LOG_DEBUG("p[%2d] = %d\n", i, p3[st[0]][st[1]][st[2] + i]);

  /* Write 3D (view) */
  if (sm_rank == 0)
    for (int i = st2[0]; i <= en2[0]; i++)
      for (int j = st2[1]; j <= en2[1]; j++)
          for (int k = st2[2]; k <= en2[2]; k++)
            q3[i][j][k] = 100 * i + 10 * j + k;

  /* Read 1D */
  for (int i = 0; i < n0 * n1 * n2; i++)
    LOG_DEBUG("q[%2d] = %d\n", i, p3[st[0]][st[1]][st[2] + i]);
  freeTensor3Shared(p3, st, sizeof(int));
}
