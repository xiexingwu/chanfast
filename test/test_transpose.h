#include <check.h>
#include <mpi_plan.h>
#include <decomp.h>
#include <transpose.h>
#include <tensor.h>
#include <io.h>

void testZY()
{
  decomp_t *decomp = createDecomp(cart_dims, cart_coord, nx, ny, nz);
  pencil_t *pz = decomp->pz;
  pencil_t *py = decomp->py;
  pencil_t *px = decomp->px;
  pencil_t *pZ = decomp->pZ;
  LOG_DEBUG("size of int: %d\n", sizeof(int));

  LOG_DEBUG("z: %d-%d, %d-%d, %d-%d\n", pz->st[0], pz->en[0], pz->st[1], pz->en[1], pz->st[2], pz->en[2]);
  LOG_DEBUG("y: %d-%d, %d-%d, %d-%d\n", py->st[0], py->en[0], py->st[1], py->en[1], py->st[2], py->en[2]);

  LOG_DEBUG("Alloc p3z\n");
  int ***p3z = tensor3(pz->st, pz->en, sizeof(int));
  LOG_DEBUG("Alloc q3y\n");
  int ***q3y = tensor3(py->st, py->en, sizeof(int));

  int *p = &p3z[pz->st[0]][pz->st[1]][pz->st[2]];
  int *q = &q3y[py->st[0]][py->st[1]][py->st[2]];

  LOG_DEBUG("Get Views\n");
  int ***p3y = tensor3View(
      &p3z[pz->st[0]][pz->st[1]][pz->st[2]],
      py->st, py->en, sizeof(int));
  int ***q3z = tensor3View(
      &q3y[py->st[0]][py->st[1]][py->st[2]],
      pz->st, pz->en, sizeof(int));

  transpose_t *trans = createTranspose(MPI_INT, decomp,
                                       &p3z[pz->st[0]][pz->st[1]][pz->st[2]],
                                       &q3y[py->st[0]][py->st[1]][py->st[2]],
                                       NULL, NULL);

  /* Write p3z */
  // for (int i = 0; i < pz->sz[0] * pz->sz[1] * pz->sz[2]; i++)
  //   p[i] = i;
  for (int i = pz->st[0]; i <= pz->en[0]; i++)
    for (int j = pz->st[1]; j <= pz->en[1]; j++)
      for (int k = pz->st[2]; k <= pz->en[2]; k++)
        p3z[i][j][k] = 100 * i + 10 * j + k;

  /* DEBUG Read */
  // for (int col = 0; col < pcol; col++)
  // {
  //   LOG_DEBUG("Send to %d\n", sm_cart_ranks[col]);
  //   for (int i = pz->st[0]; i <= pz->en[0]; i++)
  //     for (int j = pz->st[1]; j <= pz->en[1]; j++)
  //       for (int k = py->st1[col]; k <= py->en1[col]; k++)
  //         LOG_DEBUG("p3z[%d][%d][%d] = %d\n", i, j, k, p3z[i][j][k]);
  // }
  ioWrite(p3z, "p3z.dat", pz);

  MPI_Barrier(MPI_COMM_WORLD);
  LOG_DEBUG("Transpose zy\n");
  /* Transpose to q3y */
#ifdef MPI4
  zyTranspose(trans);
#else
  zyTranspose(&p3z[pz->st[0]][pz->st[1]][pz->st[2]],
              &q3y[py->st[0]][py->st[1]][py->st[2]],
              trans);
#endif

  LOG_DEBUG("READ q3y\n");
  /* Read */
  // for (int i = 0; i < py->sz[0] * py->sz[1] * py->sz[2]; i++)
  //   LOG_INFO("q[%d] = %d\n", i, q[i]);

  ioWriteY(q3y, "q3y.dat", py);
  for (int col = 0; col < pcol; col++)
  {
    LOG_DEBUG("Recv from %d\n", sm_cart_ranks[col]);
    LOG_DEBUG("%d-%d %d-%d %d-%d\n", py->st[0], py->en[0], pz->st1[col], pz->en1[col], py->st[2], py->en[2]);
    for (int i = py->st[0]; i <= py->en[0]; i++)
      for (int j = pz->st1[col]; j <= pz->en1[col]; j++)
        for (int k = py->st[2]; k <= py->en[2]; k++)
          LOG_INFO("q3y[%d][%d][%d] = %d\n", i, j, k, q3y[i][j][k]);
  }


  LOG_DEBUG("Writing q3y\n");
  /* Write q3y */
  for (int i = py->st[0]; i <= py->en[0]; i++)
    for (int j = py->st[1]; j <= py->en[1]; j++)
      for (int k = py->st[2]; k <= py->en[2]; k++)
        q3y[i][j][k] = -(100 * i + 10 * j + k);

  /* DEBUG Read */
  for (int i = pz->st[0]; i <= pz->en[0]; i++)
    for (int j = pz->st[1]; j <= pz->en[1]; j++)
      for (int k = pz->st[2]; k <= pz->en[2]; k++)
        LOG_DEBUG("q3z[%d][%d][%d] = %d\n", i, j, k, q3z[i][j][k]);

  LOG_DEBUG("Transpose yz\n");
  /* Transpose to p3z */
#ifdef MPI4
  yzTranspose(trans);
#else
  yzTranspose(&q3y[py->st[0]][py->st[1]][py->st[2]],
              &p3z[pz->st[0]][pz->st[1]][pz->st[2]],
              trans);
#endif

  /* Read */
  for (int i = pz->st[0]; i <= pz->en[0]; i++)
    for (int j = pz->st[1]; j <= pz->en[1]; j++)
      for (int k = pz->st[2]; k <= pz->en[2]; k++)
        LOG_INFO("p3z[%d][%d][%d] = %d\n", i, j, k, p3z[i][j][k]);

  freeTensor3View(p3y, py->st, sizeof(int));
  freeTensor3View(q3z, pz->st, sizeof(int));

  freeTensor3(p3z, pz->st, sizeof(int));
  freeTensor3(q3y, py->st, sizeof(int));

  MPI_Barrier(MPI_COMM_WORLD);

  freeTranspose(trans);
  freeDecomp(decomp);
}

void testZYSM()
{
  decomp_t *decomp = createDecomp(cart_dims, cart_coord, nx, ny, nz);
  transpose_t *trans = createTransposeSM(MPI_INT, decomp);
}