#include <check.h>
#include <unistd.h>


#define NARGS 2 /* Number of positional arguments */

struct args
{
  int args[NARGS]; /* it, nt */
  int prow, pcol;
  char *outputdir;
};

void helpMessage(){
  fprintf(stderr, "Usage:\n"
    "\t./a.out it nt [-r prow] [-c pcol] [-o outputdir]\n"
  );
}

int parseArgs(int argc, char **argv, struct args *args)
{
  /* Default */
  args->prow = -1;
  args->pcol = -1;
  args->outputdir = "outputdir";

  /* Positional args */
  if (argc <= 1+NARGS) {
    fprintf(stderr, "Not enough input arguments.\n");
    helpMessage();
    return 1;
  }

  for (int i = 1; i <= NARGS; i++){
    args->args[i-1] = atoi(argv[i]);
  }

  /* Optional args */
  argc -= NARGS;
  argv += NARGS;
  if (argc > 1 && argv[1][0] != '-') {
    fprintf(stderr, "Too many input arguments.\n");
    helpMessage();
    return 1;
  }

  DEBUG_PRINT0("argc: %d, argv[0]:%s\n", argc, argv[0]);
  int c;
  while ((c = getopt(argc, argv, "r:c:o:h")) != -1) {
    switch (c) {
      case 'r':
        args->prow = atoi(optarg);
        break;
      case 'c':
        args->pcol = atoi(optarg);
        break;
      case 'o':
        args->outputdir = optarg;
        break;
      case 'h':
      default:
        helpMessage();
        return 1;
    }
  }


  return 0;
}

