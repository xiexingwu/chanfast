#include <check.h>
#include <logging.h>
#include <decomp.h>
#include <globals_sim.h>
#include <io.h>

#define NULLFILE -999
#define BUFSZ 20000 // Approx 1300 numbers formatted by %+e 
#define MAXFILES 32
#define MAXLINES 100

static void appendBuf(double *u, int n, char *dst);
static void bufToAsc(char *src, char *filename);
static int fileFile(char *filename);

static char *fnames[MAXFILES];
static char *buf[MAXFILES];
static int io_round[MAXFILES]; // Number of entries in buf.

void ioAppendAscii(double u[], int n, char filename[])
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {

        int file_no = fileFile(filename); // determine which file we're opening

        appendBuf(u, n, buf[file_no]);
        io_round[file_no]+=1;

        if (io_round[file_no] == MAXLINES){
            bufToAsc(buf[file_no], fnames[file_no]);
            io_round[file_no] = 0;
            buf[file_no][0] = '\0'; // reset buffer
        }

    }
}


void ioWriteAscii(double u[], int n, char filename[])
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        FILE *fp = fopen(filename, "w");
        check(fp != NULL);
        int i;
        for (i = 0; i < n; i++) fprintf(fp, "%+e\n", u[i]);
        fclose(fp);
    }
}

void ioInitAscii(){
    LOG_WARN(0, "IO init should commit filetype/memtype for parallel io.\n");
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        int i;
        for (i = 0; i < MAXFILES; i++){
            io_round[i] = NULLFILE;
        }
    }
}
void ioFreeAscii(){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        int i;
        for (i = 0; i < MAXFILES; i++){
            if (io_round[i] != NULLFILE){
                free(fnames[i]);
                free(buf[i]);
            }
        }
    }
}
void ioFlushAscii(){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        int i;
        for (i = 0; i < MAXFILES; i++){
            if (io_round[i] == NULLFILE){ // All files done when we encounter uninitialised io_round
                break;
            }else if (io_round[i] != 0){ // non-empty buffer
                bufToAsc(buf[i], fnames[i]);
                io_round[i] = 0;
                buf[i][0] = '\0'; // reset buffer
            }
        }
    }

}

void appendBuf(double *u, int n, char *dst){
    int i;
    dst += strlen(dst);
    for (i = 0; i < n; i++){
        dst += sprintf(dst, "%+e ", u[i]);
    }
    sprintf(dst, "\n");
}
void bufToAsc(char *src, char *filename){
    FILE *fp = fopen(filename, "a");
    check(fp != NULL);
    fprintf(fp, "%s", src);
    fclose(fp);
}
int fileFile(char *filename){
    int i;
    for (i = 0; i < MAXFILES; i++){
        if (io_round[i] == NULLFILE){ // Found available index for new file
            fnames[i] = malloc(256 * sizeof(char)); // Filename should be less than 256 chars
            buf[i] = malloc(BUFSZ * sizeof(char));
            buf[i][0] = '\0';
            strcpy(fnames[i], filename);
            io_round[i] = 0;
            return i;
        }else if (strcmp(filename, fnames[i]) == 0){ // Found existing file
            return i;
        }
    }
    printf("Too many files to buffer, maximum is %d\n", MAXFILES);
    MPI_Finalize();
    exit(1);
}


void ioWrite(double ***u, char filename[],
              pencil_t *pz)
{
    LOG_WARN(0, "ioWrite: no nghostz offset in memtype.\n");
    LOG_WARN(0, "ioWrite: INT NOT DOUBLE.\n");
    MPI_Datatype filetype;
    MPI_Type_create_subarray(3,
                             (int[]){nx, ny, nz},
                             (int[]){pz->sz[0], pz->sz[1], nz},
                             (int[]){pz->st[0], pz->st[1], 0},
                             MPI_ORDER_C, MPI_INT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_Datatype memtype;
    MPI_Type_create_subarray(3,
                             (int[]){pz->sz[0], pz->sz[1], pz->sz[2]},
                             (int[]){pz->sz[0], pz->sz[1], nz},
                             (int[]){0, 0, 0},
                            //  (int[]){0, 0, nghost_z},
                             MPI_ORDER_C, MPI_INT, &memtype);
    MPI_Type_commit(&memtype);

    MPI_File fh;

    // Delete previous file if it exists
    check(MPI_SUCCESS == MPI_File_open(MPI_COMM_WORLD, filename,
        MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_DELETE_ON_CLOSE,
        MPI_INFO_NULL, &fh));

    check(MPI_SUCCESS == MPI_File_close(&fh));

    // Now write file
    check(MPI_SUCCESS == MPI_File_open(MPI_COMM_WORLD, filename,
        MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, 
        MPI_INFO_NULL, &fh));

    check(MPI_SUCCESS == MPI_File_set_view(fh, 0, MPI_INT, filetype,
        "native", MPI_INFO_NULL));

    MPI_Status status;

    check(MPI_SUCCESS == MPI_File_write_all(fh, &(u[pz->st[0]][pz->st[1]][pz->st[2]]),
        1, memtype, &status));

    check(MPI_SUCCESS == MPI_File_close(&fh));

    MPI_Type_free(&filetype);
    MPI_Type_free(&memtype);
}

void ioRead(double ***u, char filename[],
            pencil_t *pz)
{
    MPI_Datatype filetype;
    MPI_Type_create_subarray(3,
        (int []) {nx, ny, nz},
        (int []) {pz->sz[0], pz->sz[1], nz},
        (int []) {pz->st[0] - 1, pz->st[1] - 1, 0},
        MPI_ORDER_C, MPI_INT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_Datatype memtype;
    MPI_Type_create_subarray(3,
        (int []) {pz->sz[0], pz->sz[1], pz->sz[2]},
        (int []) {pz->sz[0], pz->sz[1], nz},
        (int []) {0, 0, nghost_z},
        MPI_ORDER_C, MPI_INT, &memtype);
    MPI_Type_commit(&memtype);

  //  MPI_Type_create_subarray(3,
  //      (int []) {zsz[0], zsz[1], zsz[2]},
  //      (int []) {zsz[0], zsz[1], nz},
  //      (int []) {0, 0, 0},
  //      MPI_ORDER_C, MPI_INT, &memtype);
  //  MPI_Type_commit(&memtype);

    MPI_File fh;

    check(MPI_SUCCESS == MPI_File_open(MPI_COMM_WORLD, filename,
        MPI_MODE_RDONLY, MPI_INFO_NULL, &fh));

    check(MPI_SUCCESS == MPI_File_set_view(fh, 0, MPI_INT, filetype,
        "native", MPI_INFO_NULL));

    MPI_Status status;

    check(MPI_SUCCESS == MPI_File_read_all(fh, &u[pz->st[0]][pz->st[1]][pz->st[2]],
        1, memtype, &status));

    check(MPI_SUCCESS == MPI_File_close(&fh));

    MPI_Type_free(&filetype);
    MPI_Type_free(&memtype);
}

/* DEBUG */
void ioWriteY(double ***u, char filename[],
             pencil_t *py)
{
    MPI_Datatype filetype;
    MPI_Type_create_subarray(3,
                             (int[]){nx, ny, nz},
                             (int[]){py->sz[0], py->sz[1], py->sz[2]},
                             (int[]){py->st[0], 0, py->st[2]},
                             MPI_ORDER_C, MPI_INT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_Datatype memtype;
    MPI_Type_create_subarray(3,
                             (int[]){py->sz[0], py->sz[1], py->sz[2]},
                             (int[]){py->sz[0], py->sz[1], py->sz[2]},
                             (int[]){0, 0, 0},
                             //  (int[]){0, 0, nghost_z},
                             MPI_ORDER_C, MPI_INT, &memtype);
    MPI_Type_commit(&memtype);

    MPI_File fh;

    // Delete previous file if it exists
    check(MPI_SUCCESS == MPI_File_open(MPI_COMM_WORLD, filename,
                                       MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_DELETE_ON_CLOSE,
                                       MPI_INFO_NULL, &fh));

    check(MPI_SUCCESS == MPI_File_close(&fh));
    // Now write file
    check(MPI_SUCCESS == MPI_File_open(MPI_COMM_WORLD, filename,
                                       MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL,
                                       MPI_INFO_NULL, &fh));

    check(MPI_SUCCESS == MPI_File_set_view(fh, 0, MPI_INT, filetype,
                                           "native", MPI_INFO_NULL));

    MPI_Status status;

    check(MPI_SUCCESS == MPI_File_write_all(fh, &(u[py->st[0]][py->st[1]][py->st[2]]),
                                            1, memtype, &status));

    check(MPI_SUCCESS == MPI_File_close(&fh));

    MPI_Type_free(&filetype);
    MPI_Type_free(&memtype);
}