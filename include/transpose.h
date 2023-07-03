#include <decomp.h>

/**
 * With MPI_<collective>_init, no need to expose transpose_t, just init/finalise.
 * However, will need buffers to provided for init
*/

typedef struct transpose transpose_t;

extern transpose_t transp_d; // double
extern transpose_t transp_c; // complex (spectra)

transpose_t *createTranspose(MPI_Datatype dtype, decomp_t *d,
                             void *zdata, void *ydata, void *xdata, void *Zdata);
transpose_t *createTransposeSM(MPI_Datatype type, decomp_t *decomp);
void freeTranspose(transpose_t *trans);

#ifdef MPI4
void zyTranspose(transpose_t *t);
void yzTranspose(transpose_t *t);
#else
void zyTranspose(void *src, void *dst, transpose_t *t);
void yzTranspose(void *src, void *dst, transpose_t *t);
#endif

void zyTransposeSM(void *src, void *dst, transpose_t *t);

void initTransposePlan();
void freeTransposePlan();
