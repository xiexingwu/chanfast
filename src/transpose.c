#include <check.h>
#include <logging.h>
#include <decomp.h>
#include <transpose.h>
#include <mpi_plan.h>
#include <globals_sim.h>

#ifdef MPI3
#define ALLTOALL MPI_Neighbor_alltoallw
#define ALLTOALL_init MPI_Neighbor_alltoallw_init
// #define ALLTOALL MPI_Alltoallw
// #define ALLTOALL_init MPI_Alltoallw_init
#else
#define ALLTOALL MPI_Alltoallw
#define ALLTOALL_init MPI_Alltoallw_init
#endif

/** For maximised shared memory usage (assuming contig in dim[1]/pcol), CHECK:
 * z-pencil -> x-pencil : +90 rot about y
 * x-pencil -> y-pencil : -90 rot about z
 * y-pencil -> Z-pencil : -90 rot about x
 */

typedef struct prep
{
  int *cnts;
  int *disps;
  MPI_Datatype *types0, *types1;
  MPI_Comm comm;
  MPI_Request req0, req1;
} prep_t;

typedef struct transpose
{
  prep_t zy, yx, xZ;
} transpose_t;

// static prepTranspose(prep_t *p, pencil_t *p0, pencil_t *p1, 
// MPI_Comm comm, int dim, MPI_Datatype type, void *data0, void *data1)
// {
//   /* transpose across procs with diff coord[1-dim] */
//   p->comm = comm;

//   int sz;
//   MPI_Comm_size(p->comm, &sz);
//   for (int rank = 0; rank < sz; rank++)
//   {
//     int coord[2];
//     MPI_Cart_coords(cart_comm, rank, 2, coord);
//     int rc = coord[dim];

//     p->disps[rank] = 0; /* No need for disp using subarray */

//     /* proc on different row/col of procs*/
//     if (rc != cart_coord[dim])
//     {
//       p->cnts[rank] = 0;
//       LOG_DEBUG("NULL subarray with %d.\n", rank);
//       MPI_Type_create_subarray(3,
//                                (int[]){1, 1, 1},
//                                (int[]){1, 1, 1},
//                                (int[]){0, 0, 0},
//                                MPI_ORDER_C, type, &p->types0[rank]);
//       MPI_Type_create_subarray(3,
//                                (int[]){1, 1, 1},
//                                (int[]){1, 1, 1},
//                                (int[]){0, 0, 0},
//                                MPI_ORDER_C, type, &p->types1[rank]);
//     }
//     else
//     {
//       p->cnts[rank] = 1;
//       LOG_DEBUG("zy transpose with %d.\n", rank);

//       int zsubsz[3] = {z->sz[0], z->sz[1], y->dist1[k]};
//       int zoffset[3] = {0, 0, y->st1[k]};
//       LOG_DEBUG("{%d,%d,%d}{%d,%d,%d} zst=%d\n", z->sz[0], z->sz[1], z->sz[2], zsubsz[0], zsubsz[1], zsubsz[2], zoffset[2]);
//       MPI_Type_create_subarray(3, z->sz, zsubsz, zoffset,
//                                MPI_ORDER_C, type, &zy->types0[rank]);

//       int ysubsz[3] = {y->sz[0], z->dist1[k], y->sz[2]};
//       int yoffset[3] = {0, z->st1[k], 0};
//       LOG_DEBUG("{%d,%d,%d}{%d,%d,%d} yst=%d\n", y->sz[0], y->sz[1], y->sz[2], ysubsz[0], ysubsz[1], ysubsz[2], yoffset[1]);
//       MPI_Type_create_subarray(3, y->sz, ysubsz, yoffset,
//                                MPI_ORDER_C, type, &zy->types1[rank]);
//     }

//     LOG_DEBUG("Commit subarray.\n");
//     MPI_Type_commit(&zy->types0[rank]);
//     MPI_Type_commit(&zy->types1[rank]);
//   }
// #ifdef MPI4
//   ALLTOALL_init(zdata, zy->cnts, zy->disps, zy->types0,
//                 ydata, zy->cnts, zy->disps, zy->types1,
//                 zy->comm, MPI_INFO_NULL, &zy->req0);
//   ALLTOALL_init(ydata, zy->cnts, zy->disps, zy->types1,
//                 zdata, zy->cnts, zy->disps, zy->types0,
//                 zy->comm, MPI_INFO_NULL, &zy->req1);
// #endif
// }


/** z <-> y
 * (x0, y1) <-> (x0, z1)
 */
static void zyPrep(prep_t *zy, pencil_t *z, pencil_t *y,
                   MPI_Datatype type, void *zdata, void *ydata)
{
  zy->comm = cart_comm;
  LOG_WARN(0, "zyPrep assumes x-range for all procs with same dim[0] are same.\n");

  for (int i = 0; i < cart_dims[0]; i++)
    for (int k = 0; k < cart_dims[1]; k++)
    {
      int rank;
      MPI_Cart_rank(cart_comm, (int[]){i, k}, &rank);
      zy->disps[rank] = 0; /* No need for disp using subarray */

      /* proc on different row of procs*/
      if (i != cart_coord[0])
      {
        zy->cnts[rank] = 0;
        LOG_DEBUG("NULL subarray with %d.\n", rank);
        MPI_Type_create_subarray(3,
                                 (int[]){1, 1, 1},
                                 (int[]){1, 1, 1},
                                 (int[]){0, 0, 0},
                                 MPI_ORDER_C, type, &zy->types0[rank]);
        MPI_Type_create_subarray(3,
                                 (int[]){1, 1, 1},
                                 (int[]){1, 1, 1},
                                 (int[]){0, 0, 0},
                                 MPI_ORDER_C, type, &zy->types1[rank]);
      }
      else
      {
        zy->cnts[rank] = 1;
        LOG_DEBUG("zy transpose with %d.\n", rank);

        int zsubsz[3] = {z->sz[0], z->sz[1], y->dist1[k]};
        int zoffset[3] = {0, 0, y->st1[k]};
        LOG_DEBUG("{%d,%d,%d}{%d,%d,%d} zst=%d\n", z->sz[0], z->sz[1], z->sz[2], zsubsz[0], zsubsz[1], zsubsz[2], zoffset[2]);
        MPI_Type_create_subarray(3, z->sz, zsubsz, zoffset,
                                 MPI_ORDER_C, type, &zy->types0[rank]);

        int ysubsz[3] = {y->sz[0], z->dist1[k], y->sz[2]};
        int yoffset[3] = {0, z->st1[k], 0};
        LOG_DEBUG("{%d,%d,%d}{%d,%d,%d} yst=%d\n", y->sz[0], y->sz[1], y->sz[2], ysubsz[0], ysubsz[1], ysubsz[2], yoffset[1]);
        MPI_Type_create_subarray(3, y->sz, ysubsz, yoffset,
                                 MPI_ORDER_C, type, &zy->types1[rank]);
      }

      LOG_DEBUG("Commit subarray.\n");
      MPI_Type_commit(&zy->types0[rank]);
      MPI_Type_commit(&zy->types1[rank]);
    }
#ifdef MPI4
  ALLTOALL_init(zdata, zy->cnts, zy->disps, zy->types0,
                ydata, zy->cnts, zy->disps, zy->types1,
                zy->comm, MPI_INFO_NULL, &zy->req0);
  ALLTOALL_init(ydata, zy->cnts, zy->disps, zy->types1,
                zdata, zy->cnts, zy->disps, zy->types0,
                zy->comm, MPI_INFO_NULL, &zy->req1);
#endif
}

static void zyPrepNeighbor(prep_t *zy, pencil_t *z, pencil_t *y,
                           MPI_Datatype type, void *zdata, void *ydata)
{
  zy->comm = prow_comm;
  LOG_WARN(0, "zyPrepNeighbor assumes x-range for all procs in prow_comm.\n");

  for (int rank = 0; rank < pcol; rank++)
  {
    int *coord = sm_cart_coords[rank];
    // int coord[2];
    // MPI_Cart_coords(cart_comm, sm_cart_ranks[rank], 2, coord);
    check(coord[0] == cart_coord[0]);
    int col = coord[1];

    zy->cnts[rank] = 1;
    LOG_DEBUG("zy transpose with col %d.\n", col);

    int zsubsz[3] = {z->sz[0], z->sz[1], y->dist1[col]};
    int zoffset[3] = {0, 0, y->st1[col]};
    LOG_DEBUG("{%d,%d,%d}{%d,%d,%d} zst=%d\n", z->sz[0], z->sz[1], z->sz[2], zsubsz[0], zsubsz[1], zsubsz[2], zoffset[2]);
    MPI_Type_create_subarray(3, z->sz, zsubsz, zoffset,
                             MPI_ORDER_C, type, &zy->types0[rank]);

    int ysubsz[3] = {y->sz[0], z->dist1[col], y->sz[2]};
    int yoffset[3] = {0, z->st1[col], 0};
    LOG_DEBUG("{%d,%d,%d}{%d,%d,%d} yst=%d\n", y->sz[0], y->sz[1], y->sz[2], ysubsz[0], ysubsz[1], ysubsz[2], yoffset[1]);
    MPI_Type_create_subarray(3, y->sz, ysubsz, yoffset,
                             MPI_ORDER_C, type, &zy->types1[rank]);

    LOG_DEBUG("Commit subarray.\n");
    MPI_Type_commit(&zy->types0[rank]);
    MPI_Type_commit(&zy->types1[rank]);
  }
#ifdef MPI4
  ALLTOALL_init(zdata, zy->cnts, zy->disps, zy->types0,
                ydata, zy->cnts, zy->disps, zy->types1,
                zy->comm, MPI_INFO_NULL, &zy->req0);
  ALLTOALL_init(ydata, zy->cnts, zy->disps, zy->types1,
                zdata, zy->cnts, zy->disps, zy->types0,
                zy->comm, MPI_INFO_NULL, &zy->req1);
#endif
}

static void initPrep(prep_t *p)
{
  p->cnts = malloc(world_size * sizeof(int));
  p->disps = malloc(world_size * sizeof(int));
  p->types0 = malloc(world_size * sizeof(MPI_Datatype));
  p->types1 = malloc(world_size * sizeof(MPI_Datatype));
}

static void freePrep(prep_t *p)
{
  free(p->cnts);
  free(p->disps);

  int size;
  MPI_Comm_size(p->comm, &size);
  for (int i = 0; i < size; i++)
  {
    MPI_Type_free(&p->types0[i]);
    MPI_Type_free(&p->types1[i]);
  }
  free(p->types0);
  free(p->types1);

#ifdef MPI4
  MPI_Request_free(&p->req0);
  MPI_Request_free(&p->req1);
#endif
}

transpose_t *createTranspose(MPI_Datatype dtype, decomp_t *d,
                             void *zdata, void *ydata, void *xdata, void *Zdata)
{
  transpose_t *t = malloc(sizeof(transpose_t));
  check(t != NULL);

  LOG_DEBUG("Pencils...\n");

  /** For each of x,y,z,Z pencil, need to get
   * count & displacements in row/col
   */
  initPrep(&t->zy);
  initPrep(&t->yx);
  initPrep(&t->xZ);

  LOG_DEBUG("start zyPrep\n");
#ifdef MPI3
  zyPrepNeighbor(&t->zy, d->pz, d->py, MPI_DOUBLE, zdata, ydata);
  // yxPrepNeighbor(&t->yx, d->py, d->px, MPI_DOUBLE, ydata, xdata);
  // xZPrepNeighbor(&t->xZ, d->px, d->pZ, MPI_DOUBLE, xdata, Zdata);
#else
  zyPrep(&t->zy, d->pz, d->py, MPI_DOUBLE, zdata, ydata);
  // yxPrep(&t->yx, d->py, d->px, MPI_DOUBLE, ydata, xdata);
  // xZPrep(&t->xZ, d->px, d->pZ, MPI_DOUBLE, xdata, Zdata);
#endif
  LOG_DEBUG("Finish zyPrep\n");

  return t;
}

transpose_t *createTransposeNeighbor(MPI_Datatype type, decomp_t *decomp)
{
  transpose_t *t = malloc(sizeof(transpose_t));
  check(t != NULL);
  return t;
}
transpose_t *createTransposeSM(MPI_Datatype type, decomp_t *decomp)
{
  transpose_t *t = malloc(sizeof(transpose_t));
  check(t != NULL);
  return t;
}

void freeTranspose(transpose_t *t)
{
  LOG_DEBUG("Free transpose\n");

  freePrep(&t->zy);
  // freePrep(&t->yx);
  // freePrep(&t->xZ);
  free(t);
}

#ifdef MPI4

void zyTranspose(transpose_t *t)
{
  MPI_Start(&t->zy.req0);
  MPI_Wait(&t->zy.req0, MPI_STATUS_IGNORE);
}
void yzTranspose(transpose_t *t)
{
  MPI_Start(&t->zy.req1);
  MPI_Wait(&t->zy.req1, MPI_STATUS_IGNORE);
}

#else

void zyTranspose(void *src, void *dst, transpose_t *t)
{
  prep_t *p = &t->zy;
  LOG_DEBUG("%p <-> %p\n", src, dst);
  ALLTOALL(src, p->cnts, p->disps, p->types0,
           dst, p->cnts, p->disps, p->types1, p->comm);
}
void yzTranspose(void *src, void *dst, transpose_t *t)
{
  prep_t *p = &t->zy;
  LOG_DEBUG("%p <-> %p\n", src, dst);
  ALLTOALL(src, p->cnts, p->disps, p->types1,
           dst, p->cnts, p->disps, p->types0, p->comm);
}

#endif

// void zyTransposeSM(void *src, void *dst, transpose_t *t)
// {
//   /* sm_leads with same coord[0] */
//   MPI_Neighbor_Alltoallw(src, t->zcnts_zy, t->zdispls_zy, t->ztypes_zy,
//                          dst, t->ycnts_zy, t->ydispls_zy, t->ytypes_zy, bridge_comm);
// }

void initTransposePlan()
{
  MPI_Datatype mpi_fftw_complex;
  MPI_Type_contiguous(2, MPI_DOUBLE, &mpi_fftw_complex);
  MPI_Type_commit(&mpi_fftw_complex);

  transpose_t *transp_d = createTranspose(MPI_DOUBLE, decomp_d, work1, work2, work1, work2);
  transpose_t *transp_c = createTranspose(MPI_DOUBLE, decomp_d, work1, work2, work1, work2);
}
void freeTransposePlan()
{
  freeTranspose(transp_d);
  freeTranspose(transp_c);
}
