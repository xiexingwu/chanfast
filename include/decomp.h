
#ifndef DECOMP_H
#define DECOMP_H

typedef struct pencil
{
    /* Local st/en/sz */
    int st[3], en[3], sz[3];

    /* st/en/sz of procs along row & col*/
    int *st0, *st1;
    int *en0, *en1;
    int *sz0, *sz1;
    int *dist0, *dist1;
} pencil_t;

typedef struct decomp {
    int *dims;

    pencil_t *px, *py, *pz, *pZ;

} decomp_t;

extern decomp_t *decomp_d; // double
extern decomp_t *decomp_c; // complex (spectral)

decomp_t *createDecomp(int dims[2], int coord[2],
    int nx, int ny, int nz);
void freeDecomp(decomp_t *decomp);

void initDecompPlan();
void freeDecompPlan();

#endif /* DECOMP_H */