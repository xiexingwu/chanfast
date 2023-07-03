
// Pencil provides a convenient grouping of indices related to a pencil
// TODO - should pencil be fully exposed? PIMPL to hide details? 
class Pencil
{
public:
  Pencil(); // Dummy
  Pencil(
      int *sts0, int *ens0, int *szs0,
      int *sts1, int *ens1, int *szs1);
  ~Pencil();

public:
  // Local st/en/sz 
  int st[3], en[3], sz[3];

  // st/en/sz of procs along row & col
  int *sts0, *sts1;
  int *ens0, *ens1;
  int *szs0, *szs1;
};

class Decomp
{
public:
  Decomp(); // Dummy
  Decomp(
      int dims[2], int coord[2],
      int nx, int ny, int nz);
  ~Decomp();

public:
  Pencil px, py, pz, pZ;

private:
  int *x0st, *x0en, *x0dist;
  int *y0st, *y0en, *y0dist;
  int *z0st, *z0en, *z0dist;
  int *x1st, *x1en, *x1dist;
  int *y1st, *y1en, *y1dist;
  int *z1st, *z1en, *z1dist;
};
