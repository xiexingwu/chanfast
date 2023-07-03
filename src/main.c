#include <check.h>
#include <logging.h>

#include <mpi_plan.h>
#include <tensor.h>
#include <parser.h>
#include <globals_sim.h>

#include <test_mpi_plan.h>
#include <test_tensor.h>
#include <test_transpose.h>

const double lx = 2*M_PI;
const double ly = 1*M_PI;
const double lz = 1.0;

const int nx = 8;
const int ny = 8;
const int nz = 10;

const int nghost_x = 3;
const int nghost_y = 3;
const int nghost_z = 3;

const int cflmode = 1;
const double dt = 0.001, want_cfl = 0.5;

const int ibmmode = 0;

const int istat = 500;
const int iraw = 500;
/* prow x pcol may be overwritten by command-line args*/
int prow = 2;
int pcol = 2;

const int cfrmode = 0;
const double dpdx = 1.0, dpdy = 0;
const double want_uvolavg = 0, wang_vvolavg = 0;
const double uframe = 0, vrame = 0;

const int dcdxmode = 1, ctmmode = 0; // constant thermal mass
const double want_cvolavg = 0;

const double nu = 1.0/20;
const double betg_x = 0;
const double betg_y = 0;
const double betg_z = 0;
const double prandtl = 0.71;

int RANK = -1;

int main(int argc, char **argv)
{

  /* Parse args */
  MPI_Init(&argc, &argv);
  MPI_Comm_set_errhandler( MPI_COMM_WORLD , MPI_ERRORS_RETURN);
  MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
  initLogging();


  struct args args;
  if (parseArgs(argc, argv, &args)){
    return 1;
  }
  const int nt = args.args[1];
  const int it = args.args[0];
  prow = args.prow > 0 ? args.prow : prow;
  pcol = args.prow > 0 ? args.pcol : pcol;
  const char *outputdir = args.outputdir;

  LOG_STDOUT(0, "%d: prow: %d\npcol: %d\nout: %s\n", RANK, args.prow, args.pcol, args.outputdir);

  initMpiPlan();

  // testMallocShared();
  // testMallocShared2();
  // testTensor3();
  // testTensor3Shared();
  testZY();

  /* Finalise */
  freeMpiPlan();

  freeLogging();
  DEBUG_PRINT("%d - logging freed\n", RANK);
  MPI_Finalize();
  DEBUG_PRINT("%d - MPI Finalized\n", RANK);

  return 0;
}