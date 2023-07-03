
// Forward declaration
template <typename T>
class View3;

template <typename T>
class Tensor3
{
  Tensor3(); // Dummy
  Tensor3(int st[3], int en[3], bool shared = false);
  ~Tensor3();

  // public:
  //   View3<T> getView3(int st[3], int en[3]);

public:
  bool shared;
  MPI_Win win;

private:
  const View3<T> view;
  int id = -1;

public:
  // 1D and 3D indexing using T(ijk) or T(i,j,k) instead of T[ijk], T[i][j][k]
  T &operator()(int ijk)
  {
    view(ijk);
  }
  T &operator()(int i, int j, int k)
  {
    view(i, j, k);
  };
};

template <typename T>
class View3
{
  View3(); // Dummy

  // Build View from raw data
  template <typename U>
  View3(const U *p, int st[3], int en[3]);

  // Build View from Tensor3
  template <typename U>
  View3(const Tensor3<U> &t, int st[3], int en[3])
  {
    View3(&t(0), st[3], en[3]);
  };

  ~View3();

public:
  int st[3], en[3];
  int ncells;

public:
  // 1D and 3D indexing using T(ijk) or T(i,j,k) instead of T[ijk], T[i][j][k]
  T &operator()(int ijk)
  {
    p[ijk];
  }
  T &operator()(int i, int j, int k)
  {
    p3[i][j][k];
  }
  T *p, **p2, ***p3;
};