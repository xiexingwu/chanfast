#include <check.h>
#include <mpi_plan.h>
#include <globals_sim.h>

static MPI_Win sm_win; /* will need one win for each array */
static MPI_Win sm_win2; /* for indexing in 2d */
static int x_src, x_dst, y_src, y_dst;
static int local_n = 3;
static int sm_zst[2], sm_zen[2], sm_zsz[2];

void testMallocShared()
{
  DEBUG_PRINT0("about to test 1D share malloc.\n");

  MPI_Aint local_nx = 10;
  MPI_Aint disp_unit = sizeof(int);
  int *p;
  MPI_Win_allocate_shared(local_nx * disp_unit, disp_unit, MPI_INFO_NULL, sm_comm, &p, &sm_win);

  /* Test write */
  MPI_Win_fence(0, sm_win);
  for (int i = 0; i < local_nx; i++)
    p[i] = local_nx * RANK + i;

  /* Test read */
  MPI_Win_fence(0, sm_win);
  if (sm_rank == 1)
    for (int i = 0 - sm_rank * local_nx; i < (sm_size - sm_rank) * local_nx; i++)
      printf("[%d:%d] p[%d] = %d\n", RANK, sm_rank, i, p[i]);

  /* query rank 0 memory ptr */
  int *p0;
  MPI_Aint s0;
  int d0;
  MPI_Win_shared_query(sm_win, 0, &s0, &d0, &p0);

  if (sm_rank == 1)
    for (int i = 0; i < sm_size * local_nx; i++)
      printf("[%d:%d] p0[%d] = %d\n", RANK, sm_rank, i, p0[i]);

  MPI_Win_free(&sm_win);
};

void testMallocShared2()
{
  DEBUG_PRINT0("about to test 2D share malloc.\n");

  MPI_Cart_shift(cart_comm, 0, 1, &x_src, &x_dst);
  MPI_Cart_shift(cart_comm, 1, 1, &y_src, &y_dst);

  /* Allocate 2D array */
  local_n = 3;
  MPI_Aint disp_unit = sizeof(int);

  MPI_Barrier(MPI_COMM_WORLD);

  int *p, **p_2d;
  MPI_Win_allocate_shared(sm_rank == 0 ? sm_zsz[0] * sm_zsz[1] * disp_unit : 0,
                          disp_unit, MPI_INFO_NULL, sm_comm, &p, &sm_win);
  MPI_Win_allocate_shared(sm_rank == 0 ? sm_zsz[0] * sizeof(int *) : 0,
                          sizeof(int *), MPI_INFO_NULL, sm_comm, &p_2d, &sm_win2);

  for (int i = 0; i < local_n; i++)
    p_2d[i] = p + local_n * i;

  /* Share 2d array */
  MPI_Aint s0, s0_2d;
  int d0, d0_2d;
  int *p0, **p0_2d;
  MPI_Win_shared_query(sm_win, 0, &s0, &d0, &p0);
  LOG_DEBUG("p0 -> %p\n", p0);
  MPI_Win_shared_query(sm_win2, 0, &s0_2d, &d0_2d, &p0_2d);

  /* Get shared-memory indices */
  int min_col = pcol;
  int max_col = 0;
  for (int i = 0; i < sm_size; i++)
  {
    int *coord = sm_cart_coords[i];
    if (i > 0)
      check(coord[0] == coord[-2]);
    min_col = MIN(min_col, coord[1]);
    max_col = MAX(max_col, coord[1]);
  }
  sm_zst[0] = sm_cart_coords[0][0] * local_n;
  sm_zen[0] = (sm_cart_coords[0][0] + 1) * local_n - 1;
  sm_zst[1] = min_col * local_n;
  sm_zen[1] = (max_col + 1) * local_n - 1;

  sm_zsz[0] = sm_zen[0] - sm_zst[0] + 1;
  sm_zsz[1] = sm_zen[1] - sm_zst[1] + 1;

  /* DEBUG */
  // for (int i = 0; i < world_size; i++){
  //   if (i == RANK)
  //     DEBUG_PRINT("global: %d/%d, sm: %d/%d, coord: (%d, %d), z0:%d-%d, z1:%d-%d\n", RANK, world_size - 1, sm_rank, sm_size - 1, cart_coord[0], cart_coord[1], sm_zst[0], sm_zen[0], sm_zst[1], sm_zen[1]);
  // }

  /* Root writes */
  // if (sm_rank == 0)
  // {
  //   for (int i = 0; i < sm_zsz[0] * sm_zsz[1]; i++)
  //     p[i] = i;
  //   for (int i = 0; i < sm_zsz[0]; i++)
  //     p_2d[i] = (int *)i;
  // }

  /* Each write */
  int zst[2], zen[2];
  zst[0] = sm_zst[0];
  zst[1] = sm_zst[1] + local_n * sm_rank;
  zen[0] = sm_zen[0];
  zen[1] = sm_zen[1] - local_n * (sm_size - sm_rank - 1);
  LOG_DEBUG("[%d] - (%d-%d, %d-%d)\n", sm_rank, zst[0], zen[0], zst[1], zen[1]);
  for (int i = 0; i < sm_zsz[0]; i++)
  {
    p0_2d[i] = (int *)i;
    for (int j = zst[1]; j <= zen[1]; j++)
      p0[i * sm_zsz[1] + j] = RANK * 100 + i * 10 + j;
  }

  /* read */
  MPI_Win_fence(0, sm_win);
  // MPI_Win_fence(0, sm_win2); // This fence shouldn't be needed since 2d-tensor indexing array shouldn't need to change once assigned
  for (int i = 0; i < sm_zsz[0]; i++)
  {
    LOG_INFO("p0_2d[%d] = %p\n", i, p0_2d[i]);
    for (int j = 0; j < sm_zsz[1]; j++)
      LOG_INFO("p0[%d] = %03d\n", sm_zsz[1] * i + j, p0[sm_zsz[1] * i + j]);
  }
}
