/* Domain props */
extern const double lx, ly, lz;
extern const int nx, ny, nz;
extern const int nghost_x, nghost_y, nghost_z;

/* Simulation props */
extern const int ibmmode;
extern const int istat, iraw;
extern int prow, pcol;

/* Timestepping */
extern const int cflmode;
extern const double dt, want_cfl;

/* Forcing */
extern const int cfrmode;
extern const double dpdx, dpdy;
extern const double want_uvolavg, wang_vvolavg;
extern const double uframe, vrame;

extern const int dcdxmode, ctmmode; // constant thermal mass
extern const double want_cvolavg;

/* Fluid props */
extern const double nu;
extern const double betg_x, betg_y, betg_z;
extern const double prandtl;

/* command-line args */
extern const int it, nt;
extern char outputdir[64];
