#include <chanfast/check.h>
#include <chanfast/decomp.h>

const int D_NONE = 0, D_ROW = 1, D_COL = 2;


// Decomp
Decomp::Decomp() {} // Dummy
Decomp::Decomp(
    int dims[2], int coord[2],
    int nx, int ny, int nz)
{
  x0st = new int[dims[0]];
  x0en = new int[dims[0]];
  x0dist = new int[dims[0]];
  y0st = new int[dims[0]];
  y0en = new int[dims[0]];
  y0dist = new int[dims[0]];
  z0st = new int[dims[0]];   // Unused
  z0en = new int[dims[0]];   // Unused
  z0dist = new int[dims[0]]; // Unused
  x1st = new int[dims[1]];
  x1en = new int[dims[1]];
  x1dist = new int[dims[1]];
  y1st = new int[dims[1]];
  y1en = new int[dims[1]];
  y1dist = new int[dims[1]];
  z1st = new int[dims[1]];
  z1en = new int[dims[1]];
  z1dist = new int[dims[1]];

  //         : dist
  // z-pencil: x0, y1
  // y-pencil: x0, z1
  // x-pencil: y0, z1
  // Z-pencil: y0, x1

  distribute(nx, dims[0], x0st, x0en, x0dist);
  distribute(ny, dims[0], y0st, y0en, y0dist);
  distribute(nz, dims[0], z0st, z0en, z0dist); // Unused
  distribute(nx, dims[1], x1st, x1en, x1dist);
  distribute(ny, dims[1], y1st, y1en, y1dist);
  distribute(nz, dims[1], z1st, z1en, z1dist);

  pz = Pencil(x0st, x0en, x0dist, y1st, y1en, y1dist);
  py = Pencil(x0st, x0en, x0dist, z1st, z1en, z1dist);
  px = Pencil(y0st, y0en, y0dist, z1st, z1en, z1dist);
  pZ = Pencil(y0st, y0en, y0dist, x1st, x1en, x1dist);

  // generate partition information - starting/ending index etc.
  int pdim_z[3] = {D_ROW, D_COL, D_NONE};
  int pdim_y[3] = {D_ROW, D_NONE, D_COL};
  int pdim_x[3] = {D_NONE, D_ROW, D_COL};
  int pdim_Z[3] = {D_COL, D_ROW, D_NONE};
  partition(
      nx, ny, nz,
      pdim_z, dims, coord,
      pz.st, pz.en, pz.sz);
  partition(
      nx, ny, nz,
      pdim_y, dims, coord,
      py.st, py.en, py.sz);
  partition(
      nx, ny, nz,
      pdim_x, dims, coord,
      px.st, px.en, px.sz);
  partition(
      nx, ny, nz,
      pdim_Z, dims, coord,
      pZ.st, pZ.en, pZ.sz);
}

Decomp::~Decomp()
{
  delete[] x0st;
  delete[] x0en;
  delete[] x0dist;
  delete[] y0st;
  delete[] y0en;
  delete[] y0dist;
  delete[] z0st;   // Unused
  delete[] z0en;   // Unused
  delete[] z0dist; // Unused
  delete[] x1st;
  delete[] x1en;
  delete[] x1dist;
  delete[] y1st;
  delete[] y1en;
  delete[] y1dist;
  delete[] z1st;
  delete[] z1en;
  delete[] z1dist;
}

// Pencil
Pencil::Pencil() {} // Dummy
Pencil::Pencil(
    int *st0, int *en0, int *dist0,
    int *st1, int *en1, int *dist1)
{
  st0 = st0;
  en0 = en0;
  dist0 = dist0;
  st1 = st1;
  en1 = en1;
  dist1 = dist1;
}

Pencil::~Pencil() {}

// Distribute grid points to procs along an arbitrary dimension.
// In case of uneven dist, remainder is distributed across higher rank procs.
//
// INPUT
// @param n - number of grid points in dimension to be partitioned
// @param proc - number of procs
//
// OUTPUT
// @param st[] - starting indices of procs 0, ..., n01
// @param en[] - ending indices
// @param sz[] - sizes (redundant)
void distribute(int gsz, int proc, int st[], int en[], int sz[])
{
  // gsz = lsz * proc + extra
  int lsz = gsz / proc;
  int extra = gsz - lsz * proc;

  // lower rank procs
  st[0] = 0;
  en[0] = lsz - 1;
  sz[0] = lsz;
  for (int i = 1; i < proc - extra; i++)
  {
    st[i] = en[i - 1] + 1;
    en[i] = en[i - 1] + lsz;
    sz[i] = lsz;
  }

  // higher rank procs
  lsz = lsz + 1;
  for (int i = proc - extra; i < proc; i++)
  {
    st[i] = en[i - 1] + 1;
    en[i] = en[i - 1] + lsz;
    sz[i] = lsz;
  }
}

// Determine sub-domain info of current proc
//
// INPUT
// @param nx - global nx
// @param ny - global ny
// @param nz - global nz
// @param pdim[3] - array to indicate whether x,y,z are distributed locally or across dims[0] or dims[1]
//
// OUTPUT
// @param lstart[3] - starting indices in x,y,z
// @param lend[3] - ending indices
// @param lsz[3] - sizes (redundant)
void partition(
    int nx, int ny, int nz, int pdim[3],
    int dims[2], int coord[2],
    int lstart[3], int lend[3], int lsize[3])
{
  for (int i = 0; i < 3; i++)
  {
    int gsize;
    if (i == 0)
      gsize = nx;
    if (i == 1)
      gsize = ny;
    if (i == 2)
      gsize = nz;

    // Not distributed
    if (pdim[i] == D_NONE)
    { // all local
      lstart[i] = 0;
      lend[i] = gsize - 1;
      lsize[i] = gsize;
      return;
    }

    // Distributed
    int *st = new int[dims[0]];
    int *en = new int[dims[0]];
    int *sz = new int[dims[0]];
    if (pdim[i] == D_ROW)
    {
      // distribute across dims[0]
      distribute(gsize, dims[0], st, en, sz);
      lstart[i] = st[coord[0]];
      lend[i] = en[coord[0]];
      lsize[i] = sz[coord[0]];
    }
    else if (pdim[i] == D_COL)
    {
      // distribute across dims[1]
      distribute(gsize, dims[1], st, en, sz);
      lstart[i] = st[coord[1]];
      lend[i] = en[coord[1]];
      lsize[i] = sz[coord[1]];
    }
    delete[] st;
    delete[] en;
    delete[] sz;
  }
}
