#include <math.h>
#include <mpi.h>

#include <spdlog/spdlog.h>

#include <chanfast/params/args.h>
#include <chanfast/params/grid.h>
#include <chanfast/params/sim.h>

#include <chanfast/check.h>
#include <chanfast/tensor.h>
#include <chanfast/decomp.h>
#include <chanfast/mpi_plan.h>

namespace grid
{
  const double lx = 2 * M_PI;
  const double ly = 1 * M_PI;
  const double lz = 1.0;
  const int nx = 16;
  const int ny = 16;
  const int nz = 40;
}

namespace args
{
  const int istat = 100;
  const int iraw = 100;
  const int prow = 2;
  const int pcol = 2;
  const int it = 1;
  const int nt = 1;
  char outputdir[64] = "outputdir";
}

namespace sim
{
  const int cflmode = 1;
  const double dt = 0.001, want_cfl = 0.5;

  const int cfrmode = 0;
  const double dpdx = 1.0, dpdy = 0;
  const double want_uvolavg = 0, wang_vvolavg = 0;
  const double uframe = 0, vrame = 0;

  const int dcdxmode = 1, ctmmode = 0; // constant thermal mass
  const double want_cvolavg = 0;

  const double nu = 1.0 / 20;
  const double betg_x = 0, betg_y = 0, betg_z = 0;
  const double prandtl = 0.71;
}

int main(int argc, char **argv)
{

  /* Parse args */
  MPI_Init(&argc, &argv);

  // initLogging();

  // struct args args;
  // if (parseArgs(argc, argv, &args)){
  //   return 1;
  // }
  // const int nt = args.args[1];
  // const int it = args.args[0];
  // prow = args.prow > 0 ? args.prow : prow;
  // pcol = args.prow > 0 ? args.pcol : pcol;
  // const char *outputdir = args.outputdir;

  // DEBUG_PRINT0("%d: prow: %d\npcol: %d\nout: %s\n", RANK, args.prow, args.pcol, args.outputdir);

  Mpi mpi = Mpi();

  // // testMallocShared();
  // // testMallocShared2();
  // // testTensor3();
  // testTensor3Shared();

  // MPI_Barrier(MPI_COMM_WORLD);
  // freeLogging();
  MPI_Finalize();
  return 0;
}