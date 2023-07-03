#include <gtest/gtest.h>

#include <chanfast/decomp.h>

// TEST(DecompTest, Small)
// {
//   const int nx = 10, ny = 15, nz = 20;
//   const int prow = 3, pcol = 4;
//   int dims[2] = {prow, pcol};

//   Decomp *d = new Decomp[prow * pcol];
//   for (int i = 0; i < prow; i++)
//     for (int j = 0; j < pcol; j++)
//       d[i * prow + j] = Decomp(dims, (int[2]){i, j}, nx, ny, nz);

//   int x, y, z, X, Y, Z;

//   // z-pencil sizes
//   X = Y = Z = 0;
//   for (int i = 0; i < prow; i++)
//   {
//     printf("zsz %d %d %d\n", d[i * prow].zsz[0], d[i * prow].zsz[1], d[i * prow].zsz[2]);
//     X += d[i * prow].zsz[0];
//     for (int j = 0; j < pcol; j++)
//       EXPECT_EQ(nz, d[i * prow + j].zsz[2]);
//   }
//   for (int j = 0; j < pcol; j++)
//   {
//     printf("zsz %d %d %d\n", d[j].zsz[0], d[j].zsz[1], d[j].zsz[2]);
//     Y += d[j].zsz[1];
//   }
//   EXPECT_EQ(nx, X);
//   EXPECT_EQ(ny, Y);

//   // z-pencil indices
//   for (int i = 0; i < prow; i++)
//   {
//     EXPECT_EQ(d[i * prow].zst[0], 0);
//     EXPECT_EQ(d[i * prow + pcol - 1].zen[0], nx - 1);
//   }
//   for (int j = 0; j < pcol; j++)
//   {
//     EXPECT_EQ(d[j].zst[1], 0);
//     EXPECT_EQ(d[(pcol - 1) * prow + j].zen[1], ny - 1);
//   }

//   for (int i = 1; i < prow; i++)
//     for (int j = 1; j < pcol; j++)
//     {
//       EXPECT_EQ(d[i * prow + j].zst[2], 0);
//       EXPECT_EQ(d[i * prow + j].zst[0], d[(i - 1) * prow + j].zen[0] + 1);
//       EXPECT_EQ(d[i * prow + j].zst[1], d[i * prow + j + 1].zen[1] + 1);
//       EXPECT_EQ(d[i * prow + j].zen[2], nz - 1);
//     }
// }