#include <check.h>
#include <logging.h>
#include <decomp.h>
#include <mpi_plan.h>

#include <globals_sim.h>

enum dist_t
{
  DIST0,
  DIST1,
  DISTL
}; // flag to distribute locally, across dims[0], across dims[1]

static int *x0st, *x0en, *x0dist;
static int *y0st, *y0en, *y0dist;
static int *z0st, *z0en, *z0dist;
static int *x1st, *x1en, *x1dist;
static int *y1st, *y1en, *y1dist;
static int *z1st, *z1en, *z1dist;

static pencil_t *initPencil(int *st0, int *en0, int *dist0,
                            int *st1, int *en1, int *dist1)
{
  pencil_t *p = malloc(sizeof(pencil_t));
  p->st0 = st0;
  p->en0 = en0;
  p->dist0 = dist0;
  p->st1 = st1;
  p->en1 = en1;
  p->dist1 = dist1;
  return p;
}

/**
 * Distribute grid points to procs along an arbitrary dimension.
 * In case of uneven dist, remainder is distributed across higher rank procs.
 *
 * INPUT
 * @param n - number of grid points in dimension to be partitioned
 * @param proc - number of procs
 *
 * OUTPUT
 * @param st[] - starting indices of procs 0, ..., n01
 * @param en[] - ending indices
 * @param sz[] - sizes (redundant)
 */
static void distribute(int n, int proc, int st[], int en[], int sz[])
{
  /* n = q * nproc + r */
  int q = n / proc;
  int r = n - q * proc;

  int i;
  st[0] = 0;
  sz[0] = q;
  en[0] = q - 1;

  /* lower rank procs */
  for (i = 1; i < proc - r; i++)
  {
    st[i] = en[i - 1] + 1;
    sz[i] = q;
    en[i] = en[i - 1] + q;
  }
  /* higher rank procs */
  q = q + 1;
  for (i = proc - r; i < proc; i++)
  {
    st[i] = en[i - 1] + 1;
    sz[i] = q;
    en[i] = en[i - 1] + q;
  }
}

/**
 * Determine sub-domain info of current proc
 *
 * INPUT
 * @param nx - global nx
 * @param ny - global ny
 * @param nz - global nz
 * @param pdim[3] - array to indicate whether x,y,z are distributed locally or across dims[0] or dims[1]
 *
 * OUTPUT
 * @param lstart[3] - starting indices in x,y,z
 * @param lend[3] - ending indices
 * @param lsz[3] - sizes (redundant)
 */
static void partition(int nx, int ny, int nz, enum dist_t pdim[3],
                      int dims[2], int coord[2],
                      int lstart[3], int lend[3], int lsize[3])
{
  int *st, *en, *sz;
  int i, gsize;

  for (i = 0; i < 3; i++)
  {
    switch (i)
    {
    case 0:
      gsize = nx;
      break;
    case 1:
      gsize = ny;
      break;
    case 2:
      gsize = nz;
      break;
    }

    int d = -1; // split in dims[d]
    switch (pdim[i])
    {
    case DISTL:
      lstart[i] = 0;
      lend[i] = gsize - 1;
      lsize[i] = gsize;
      continue;
    case DIST0:
      d = 0;
      break;
    case DIST1:
      d = 1;
      break;
    }

    st = (int *)malloc(sizeof(int) * dims[d]);
    en = (int *)malloc(sizeof(int) * dims[d]);
    sz = (int *)malloc(sizeof(int) * dims[d]);
    check(st != NULL);
    check(en != NULL);
    check(sz != NULL);
    distribute(gsize, dims[d], st, en, sz);
    lstart[i] = st[coord[d]];
    lend[i] = en[coord[d]];
    lsize[i] = sz[coord[d]];
    free(st);
    free(en);
    free(sz);
  }
}

/**
 * Determine how each of nx, ny, nz are distributed across processors in each pencil
 *
 * INPUT
 * @param nx - global nx
 * @param ny - global ny
 * @param nz - global nz
 * @param dims[2] - [prow, pcol]
 *
 * OUTPUT
 * @param [xyz][01](st|en) - start end indices for n[xyz] along prow/pcol
 * @param [xyz][01]dist - nprocs  n[xyz] along prow/pcol
 */
void getDist(int nx, int ny, int nz, decomp_t *d)
{
  int *dims = d->dims;

  /**        : dist
   * z-pencil: x0, y1
   * y-pencil: x0, z1
   * x-pencil: y0, z1
   * Z-pencil: y0, x1
   */
  LOG_DEBUG("Dist nx in 0.\n");
  distribute(nx, dims[0], x0st, x0en, x0dist);
  LOG_DEBUG("Dist ny in 0.\n");
  distribute(ny, dims[0], y0st, y0en, y0dist);
  LOG_DEBUG("Dist nz in 0.\n");
  distribute(nz, dims[0], z0st, z0en, z0dist); // Unused

  LOG_DEBUG("Dist nx in 1.\n");
  distribute(nx, dims[1], x1st, x1en, x1dist);
  LOG_DEBUG("Dist ny in 1.\n");
  distribute(ny, dims[1], y1st, y1en, y1dist);
  LOG_DEBUG("Dist nz in 1.\n");
  distribute(nz, dims[1], z1st, z1en, z1dist);

  LOG_DEBUG("Finished dist.\n");
}

decomp_t *createDecomp(int dims[2], int coord[2],
                       int nx, int ny, int nz)
{
  decomp_t *d = (decomp_t *)malloc(sizeof(decomp_t));
  check(d != NULL);

  LOG_DEBUG("Starting decomp.\n");

  d->dims = cart_dims;

  // distribute mesh points
  x0st = malloc(dims[0] * sizeof(int));
  x0en = malloc(dims[0] * sizeof(int));
  x0dist = malloc(dims[0] * sizeof(int));
  y0st = malloc(dims[0] * sizeof(int));
  y0en = malloc(dims[0] * sizeof(int));
  y0dist = malloc(dims[0] * sizeof(int));
  z0st = malloc(dims[0] * sizeof(int));   // Unused
  z0en = malloc(dims[0] * sizeof(int));   // Unused
  z0dist = malloc(dims[0] * sizeof(int)); // Unused
  x1st = malloc(dims[1] * sizeof(int));
  x1en = malloc(dims[1] * sizeof(int));
  x1dist = malloc(dims[1] * sizeof(int));
  y1st = malloc(dims[1] * sizeof(int));
  y1en = malloc(dims[1] * sizeof(int));
  y1dist = malloc(dims[1] * sizeof(int));
  z1st = malloc(dims[1] * sizeof(int));
  z1en = malloc(dims[1] * sizeof(int));
  z1dist = malloc(dims[1] * sizeof(int));

  check(x0st != NULL);
  check(x0en != NULL);
  check(x0dist != NULL);
  check(y0st != NULL);
  check(y0en != NULL);
  check(y0dist != NULL);
  check(z1st != NULL);
  check(z1en != NULL);
  check(z1dist != NULL);
  check(x1st != NULL);
  check(x1en != NULL);
  check(x1dist != NULL);
  check(y1st != NULL);
  check(y1en != NULL);
  check(y1dist != NULL);
  check(z1st != NULL);
  check(z1en != NULL);
  check(z1dist != NULL);

  d->pz = initPencil(x0st, x0en, x0dist, y1st, y1en, y1dist);
  d->py = initPencil(x0st, x0en, x0dist, z1st, z1en, z1dist);
  d->px = initPencil(y0st, y0en, y0dist, z1st, z1en, z1dist);
  d->pZ = initPencil(y0st, y0en, y0dist, x1st, x1en, x1dist);

  getDist(nx, ny, nz, d);

  LOG_DEBUG("Finished dist.\n");

  // generate partition information - starting/ending index etc.
  partition(nx, ny, nz, (int[]){DIST0, DIST1, DISTL}, dims, coord,
            d->pz->st, d->pz->en, d->pz->sz);
  partition(nx, ny, nz, (int[]){DIST0, DISTL, DIST1}, dims, coord,
            d->py->st, d->py->en, d->py->sz);
  partition(nx, ny, nz, (int[]){DISTL, DIST0, DIST1}, dims, coord,
            d->px->st, d->px->en, d->px->sz);
  partition(nx, ny, nz, (int[]){DIST1, DIST0, DISTL}, dims, coord,
            d->pZ->st, d->pZ->en, d->pZ->sz);

  return d;
}

void freeDecomp(decomp_t *d)
{
  LOG_DEBUG("Freed decomp\n");
  free(d->pz);
  free(d->py);
  free(d->px);
  free(d->pZ);

  free(x0st);
  free(x0en);
  free(x0dist);
  free(y0st);
  free(y0en);
  free(y0dist);
  free(z0st);
  free(z0en);
  free(z0dist);
  free(x1st);
  free(x1en);
  free(x1dist);
  free(y1st);
  free(y1en);
  free(y1dist);
  free(z1st);
  free(z1en);
  free(z1dist);

  free(d);
}

void initDecompPlan()
{
  decomp_t *decomp_d = createDecomp(cart_dims, cart_coord, nx, ny, nz + 2 * nghost_z);
  decomp_t *decomp_c = createDecomp(cart_dims, cart_coord, nx / 2 + 1, ny, nz + 2 * nghost_z);
}

void freeDecompPlan()
{
  freeDecomp(decomp_d);
  freeDecomp(decomp_c)
}