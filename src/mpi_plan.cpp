#include <mpi.h>

#include <spdlog/spdlog.h>

#include <chanfast/check.h>
#include <chanfast/mpi_plan.h>
#include <chanfast/params/args.h>

Mpi::Mpi()
{
  world_cm = MPI_COMM_WORLD;
  MPI_Comm_size(MPI_COMM_WORLD, &world_sz);
  check(args::prow * args::pcol == world_sz);
  dims[0] = args::prow;
  dims[1] = args::pcol;
  periodic[0] = 1;
  periodic[1] = 1;

  // Cart
  MPI_Cart_create(world_cm, 2, dims, periodic, 0, &cart_cm);
  MPI_Cart_coords(cart_cm, world_rank, 2, coord);

  /* Shared-memory comm */
  MPI_Comm tmp_sm;
  MPI_Comm_split_type(world_cm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &tmp_sm);
  // Don't span multiple proc_rows
  MPI_Comm_split(tmp_sm, coord[0], coord[1], &sm_cm);
  MPI_Comm_rank(sm_cm, &sm_rank);
  MPI_Comm_size(sm_cm, &sm_sz);
  MPI_Comm_free(&tmp_sm);

  // // prow pcol sub-comms
  // MPI_Comm_split(MPI_COMM_WORLD, coord[0], coord[1], &prow_comm);
  // MPI_Comm_split(MPI_COMM_WORLD, coord[1], coord[0], &pcol_comm);
}

Mpi::~Mpi()
{
  MPI_Comm_free(&cart_cm);
}

