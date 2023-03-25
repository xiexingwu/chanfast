void ioAppendAscii(double u[], int n, char filename[]);
void ioWriteAscii(double u[], int n, char filename[]);

void ioFlushAscii();
void ioInitAscii();
void ioFreeAscii();

void ioWrite(double ***u, char filename[],
             pencil_t *pz);
void ioRead(double ***u, char filename[],
             pencil_t *pz);

/* DEBUG */
void ioWriteY(double ***u, char filename[],
              pencil_t *py);
