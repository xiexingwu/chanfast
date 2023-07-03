/* Simulation props */
namespace sim
{
  extern const int ibmmode;

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
}
