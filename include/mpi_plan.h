
extern void initMpiPlan();

extern MPI_Comm world_comm, sm_comm, cart_comm;
extern int sm_rank;
extern int world_size, sm_size;
extern int cart_dims[2], cart_coord[2], periodic[2];
extern int *sm_cart_ranks, *sm_cart_coords[2]; // cart rank/coords of all procs in sm
