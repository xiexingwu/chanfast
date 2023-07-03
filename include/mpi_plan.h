#ifndef MPI_PLAN_H
#define MPI_PLAN_H

extern void initMpiPlan();
extern void freeMpiPlan();

extern MPI_Comm world_comm, sm_comm, cart_comm, bridge_comm, prow_comm, pcol_comm;
extern int sm_rank;
extern int world_size, sm_size, bridge_size;
extern int cart_dims[2], cart_coord[2], periodic[2];
extern int *sm_cart_ranks, *sm_cart_coords[2]; // cart rank/coords of all procs in sm
extern int *bridge_ranks;

#endif /* MPI_PLAN_H */
