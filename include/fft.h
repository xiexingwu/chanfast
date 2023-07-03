// #include <decomp.h>

typedef struct fftx
{
  pencil_t *px_d;
  pencil_t *px_c;
  fftw_plan d2c;
  fftw_plan c2d;
} fftx_plan;

typedef struct ffty
{
  pencil_t *py_c;
  fftw_plan i2o;
  fftw_plan o2i;
} fftx_plan;

void initFftPlan();
void freeFftPlan();

void fftx_d2c(fftx_plan *p);
void fftx_c2d(fftx_plan *p);
void fftx_i2o(ffty_plan *p);
void fftx_o2i(ffty_plan *p);