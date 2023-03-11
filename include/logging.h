
extern void initLogging();
extern void LOG_STDOUT(int rank, const char *fmt, ...);
extern void LOG_WARN(int rank, const char *fmt, ...);
extern void LOG_ERR(int rank, const char *fmt, ...);

extern void LOG_INFO(const char *fmt, ...);
extern void LOG_DEBUG(const char *fmt, ...);
extern void freeLogging();
