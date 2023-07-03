#include <check.h>
#include <logging.h>
#include <globals_sim.h>
#include <mpi_plan.h>

/**
 * See links below for basic shared memory MPI tutorial
 * https://www.intel.com/content/dam/develop/external/us/en/documents/an-introduction-to-mpi-3-597891.pdf
 * https://fs.hlrs.de/projects/par/mooc/mooc-2/mooc2-week3-4.pdf
 */
/**
 * For MPI+MPI collectives implementations
 * https://doi.org/10.1016/j.parco.2020.102669
 */
void translateSmToCart();
void checkCartContiguous();

MPI_Comm world_comm, sm_comm, cart_comm, bridge_comm, prow_comm, pcol_comm;
int world_size, sm_size, bridge_size;
int sm_rank;
int cart_dims[2], cart_coord[2], periodic[2];
int *sm_cart_ranks, *sm_cart_coords[2]; // cart rank/coords of all procs in sm
int *bridge_ranks;

MPI_Group cart_group, sm_group;

void initMpiPlan()
{
  world_comm = MPI_COMM_WORLD;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  check(prow * pcol == world_size);
  cart_dims[0] = prow;
  cart_dims[1] = pcol;
  periodic[0] = 1;
  periodic[1] = 1;

  /* TESTING: split world in half */
  LOG_WARN(0, "Splitting MPI_COMM_WORLD...\n");
  MPI_Comm_split(MPI_COMM_WORLD, RANK < world_size / 2, RANK, &world_comm);

  /* Cartesian comm */
  LOG_WARN(0, "Cart_create assuming block-partitioning of CPUs on nodes?\n");
  MPI_Cart_create(MPI_COMM_WORLD, 2, cart_dims, periodic, 0, &cart_comm);
  MPI_Cart_coords(cart_comm, RANK, 2, cart_coord);

  /* prow pcol sub-comms */
  MPI_Comm_split(MPI_COMM_WORLD, cart_coord[0], cart_coord[1], &prow_comm);
  MPI_Comm_split(MPI_COMM_WORLD, cart_coord[1], cart_coord[0], &pcol_comm);

  /* Shared-memory comm */
  MPI_Comm tmp_sm;
  MPI_Comm_split_type(world_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &tmp_sm);
  // Don't span multiple proc_rows
  MPI_Comm_split(tmp_sm, cart_coord[0], cart_coord[1], &sm_comm);
  MPI_Comm_rank(sm_comm, &sm_rank);
  MPI_Comm_size(sm_comm, &sm_size);
  MPI_Comm_free(&tmp_sm);

  MPI_Comm_split(world_comm, sm_rank, RANK, &bridge_comm);
  MPI_Comm_size(bridge_comm, &bridge_size);
  bridge_ranks = malloc(bridge_size * sizeof(int));
  // MPI_Allgather(&RANK, 1, MPI_INT, bridge_ranks, 1, MPI_INT, bridge_comm);

  /* Check contiguity */
  translateSmToCart();
  checkCartContiguous();
}

void freeMpiPlan()
{
  LOG_DEBUG("Free MpiPlan...\n");

  free(bridge_ranks);
  free(sm_cart_ranks);
  free(&(sm_cart_coords[0][0]));  

  MPI_Comm_free(&bridge_comm);
  MPI_Comm_free(&sm_comm);

  MPI_Comm_free(&prow_comm);
  MPI_Comm_free(&pcol_comm);

  MPI_Comm_free(&cart_comm);
  MPI_Comm_free(&world_comm);
}

/**
 * Get cart_rank and cart_coords of all procs in sm_comm.
 */
void translateSmToCart()
{
  MPI_Comm_group(cart_comm, &cart_group);
  MPI_Comm_group(sm_comm, &sm_group);

  sm_cart_ranks = malloc(sm_size * sizeof(int));
  int *tmp = malloc(sm_size * 2 * sizeof(int)); // Free via sm_cart_coords
  for (int i = 0; i < sm_size; i++)
  {
    sm_cart_coords[i] = tmp + 2 * i;
    MPI_Group_translate_ranks(sm_group, 1, &i, cart_group, sm_cart_ranks + i);
    MPI_Cart_coords(cart_comm, sm_cart_ranks[i], 2, sm_cart_coords[i]);
  }
}

/**
 * Convert all sm_cart_coords to 1d indexing, and check if contiguous.
 */
void checkCartContiguous()
{
  int *inds = malloc(sm_size * sizeof(int));
  int min_ind = world_size;

  /* Get indices */
  LOG_DEBUG("Getting indices.\n");
  for (int i = 0; i < sm_size; i++)
  {
    int *coord = sm_cart_coords[i];
    if (i > 0)
      check(coord[0] == coord[-2]);
    inds[i] = (coord[0] * pcol + coord[1]) % world_size;
    min_ind = MIN(min_ind, inds[i]);
  }

  /* Check contiguous */
  LOG_DEBUG("Checking contig.\n");
  for (int ind = min_ind; ind < min_ind + sm_size; ind++)
  {
    int found = 0;
    for (int i = 0; i < sm_size; i++)
      if (inds[i] == ind)
        found = 1; // no need to break
    if (!found)
    {
      LOG_ERR(RANK, "NONCONTIG: index %d not found betwen (%d, %d).\n", ind, min_ind, min_ind + sm_size);
      abort();
    }
  }

  /* DEBUG: print coords of all procs in sm_comm */
  for (int i = 0; i < world_size; i++)
    if (i == RANK && sm_rank == 0)
      for (int j = 0; j < sm_size; j++)
        LOG_STDOUT(RANK, "\t[%d:%d] - (%d, %d)\n", sm_cart_ranks[j], j, sm_cart_coords[j][0], sm_cart_coords[j][1]);

  free(inds);
}

void initDecomp()
{
  /**
   * 1. sm_comms contiguous in p_row (arbitrary choice)
   *    i. if
   * 2. sm_comm leads need decomp_plan for the sm_comm
   *     i. add up local sizes in decomp-direction
   *    ii.
   */
  /* Decomp only for sm_comm leads */

  /* z-pencil */
  /* prow: y, pcol: x*/

  /* y-pencil */

  /* x-pencil */
  /* z-pencil2 (xz-transpose for poisson solver) */
}