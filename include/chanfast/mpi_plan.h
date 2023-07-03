
class Mpi 
{
public:
  Mpi();
  ~Mpi();

public:
  static inline MPI_Comm world_cm, sm_cm;
  static inline int world_rank, sm_rank;
  static inline int world_sz, sm_sz;

  static inline MPI_Comm cart_cm;
  static inline int dims[2], periodic[2], coord[2];
};
