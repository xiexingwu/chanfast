
void initLogging();
void LOG_STDOUT(int rank, const char *fmt, ...);
void LOG_WARN(int rank, const char *fmt, ...);
void LOG_ERR(int rank, const char *fmt, ...);
void LOG_INFO(const char *fmt, ...);
void LOG_DEBUG(const char *fmt, ...);
void freeLogging();
