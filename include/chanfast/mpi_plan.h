
class Mpi 
{
public:
  Mpi();
  ~Mpi();

public:
  static MPI_Comm world_cm, sm_cm;
  static int world_rank, sm_rank;
  static int world_sz, sm_sz;

  static MPI_Comm cart_cm;
  static int dims[2], periodic[2], coord[2];
};
