#include <stdlib.h>
#include <mpi.h>

#define MAX_WINS 64
extern MPI_Win tens_wins[MAX_WINS];
extern void *tens_adds[MAX_WINS];
extern int findAvailWin();
extern int findWinByAdd(void *p);

void ***tensor3(int st[3], int en[3], size_t dtype_sz);
void ***tensor3Shared(int st[3], int en[3], size_t dtype_sz);
void ***tensor3View(void *p, int st[3], int en[3], size_t dtype_sz);

void freeTensor3(void ***p3, int st[3], size_t dtype_sz);
void freeTensor3View(void ***p3, int st[3], size_t dtype_sz);
void freeTensor3Shared(void ***p3, int st[3], size_t dtype_sz);

void fenceTensor3(void *p);
