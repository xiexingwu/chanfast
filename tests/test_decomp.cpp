#include <gtest/gtest.h>

#include <chanfast/decomp.h>

const int nx = 10, ny = 15, nz = 20;
const int prow = 3, pcol = 4;
int dims[2] = {prow, pcol};

TEST(DecompTest, Distribute)
{
  int prow = 4;
  int nx = 15;
  int *x0st = new int[prow];
  int *x0en = new int[prow];
  int *x0sz = new int[prow];
  distribute(nx, prow, x0st, x0en, x0sz);

  EXPECT_EQ(x0sz[0], 3);
  int sz = x0sz[0];
  EXPECT_EQ(x0st[0], 0);
  EXPECT_EQ(x0st[0] + x0sz[0] - 1, x0en[0]);
  for (int i = 1; i < prow; i++)
  {
    EXPECT_EQ(x0sz[i], 4);
    sz += x0sz[i];
    EXPECT_EQ(x0st[i], x0en[i - 1] + 1);
    EXPECT_EQ(x0st[i] + x0sz[i] - 1, x0en[i]);
  }
  EXPECT_EQ(x0en[prow - 1], nx - 1);
  EXPECT_EQ(sz, nx);
}

TEST(DecompTest, Partition)
{
  int pdim[3] = {D_ROW, D_COL, D_NONE};
  int coord[2] = {prow - 1, pcol - 1};
  int lst[3], len[3], lsz[3];
  partition(nx, ny, nz, pdim, dims, coord, lst, len, lsz);
  EXPECT_EQ(lsz[0], nx / prow + (nx % prow ? 1 : 0));
  EXPECT_EQ(lsz[1], ny / pcol + (ny % pcol ? 1 : 0));
  EXPECT_EQ(lsz[2], nz);

  EXPECT_EQ(lst[0], nx - lsz[0]);
  EXPECT_EQ(lst[1], ny - lsz[1]);
  EXPECT_EQ(lst[2], 0);

  EXPECT_EQ(len[0], nx - 1);
  EXPECT_EQ(len[1], ny - 1);
  EXPECT_EQ(len[2], nz - 1);
}

TEST(DecompTest, Small)
{

  Decomp *d = new Decomp[prow * pcol];
  Pencil *pz = new Pencil[prow * pcol];
  for (int i = 0; i < prow; i++)
    for (int j = 0; j < pcol; j++)
    {
      d[i * pcol + j] = Decomp(dims, (int[2]){i, j}, nx, ny, nz);
      pz[i * pcol + j] = d[i * pcol + j].pz;
    }

  int x, y, z, X, Y, Z;

  // z-pencil sizes
  X = Y = Z = 0;
  for (int i = 0; i < prow; i++)
  {
    X += pz[i * pcol].sz[0];
    for (int j = 0; j < pcol; j++)
    {
      // printf("[%d,%d] - zsz %d %d %d\n", i, j, pz[i*pcol+j].sz[0], pz[i*pcol+j].sz[1], pz[i*pcol+j].sz[2]);
      EXPECT_EQ(nz, pz[i * pcol + j].sz[2]);
    }
  }
  for (int j = 0; j < pcol; j++)
  {
    // printf("zsz %d %d %d\n", pz[j].sz[0], pz[j].sz[1], pz[j].sz[2]);
    Y += pz[j].sz[1];
  }
  EXPECT_EQ(nx, X);
  EXPECT_EQ(ny, Y);

  // z-pencil indices
  // count along col, pz[0,j] should start at x=0, pz[prow-1,j] should end at x=nx-1
  for (int j = 0; j < pcol; j++)
  {
    EXPECT_EQ(pz[j].st[0], 0);
    EXPECT_EQ(pz[(prow - 1) * pcol + j].en[0], nx - 1);
  }
  // count along row, pz[i,0] should start at y=0, pz[i,pcol-1] should end at y=ny-1
  for (int i = 0; i < prow; i++)
  {
    EXPECT_EQ(pz[i * pcol].st[1], 0);
    EXPECT_EQ(pz[i * pcol + pcol - 1].en[1], ny - 1);
  }

  for (int i = 1; i < prow; i++)
    for (int j = 1; j < pcol; j++)
    {
      EXPECT_EQ(pz[i * pcol + j].st[0], pz[(i - 1) * pcol + j].en[0] + 1);
      EXPECT_EQ(pz[i * pcol + j].st[1], pz[i * pcol + j - 1].en[1] + 1);

      // z should be trivially 0 to nz-1
      EXPECT_EQ(pz[i * pcol + j].st[2], 0);
      EXPECT_EQ(pz[i * pcol + j].en[2], nz - 1);
    }
}