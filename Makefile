SDIR = ./src
IDIR = ./include
ODIR = ./src/obj
TDIR = ./test
LOGS = ./logs

CC = mpicc
CFLAGS = -I$(IDIR) -I$(TDIR)
LIBS = -lm -lfftw3
WARNS = -Wno-discarded-qualifiers -Wno-int-to-pointer-cast -Wno-incompatible-pointer-types

_DEP = check.h logging.h parser.h globals_sim.h fourth.h rk3.h
_OBJ = main.o logging.o mpi_plan.o tensor.o
_TST = test_mpi_plan.h test_tensor.h

TARGET = ./a.out

DEP = $(patsubst %,$(IDIR)/%,$(_DEP))
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))
TST = $(patsubst %,$(TDIR)/%,$(_TST))

$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(DEFINES) $(INCLUDES)

$(ODIR)/%.o: $(SDIR)/%.c $(DEP) $(TST)
	$(CC) -c -o $@ $< $(CFLAGS) $(DEFINES) $(INCLUDES) $(WARNS)

.PHONY: clean clobber all test
clean:
	rm -f $(ODIR)/*.o $(LOGS)/*

clobber: clean
	rm -f $(TARGET)

all: clobber $(TARGET)

test:
	make DEFINES=-DDEBUG all
	mpirun -n 4 ./a.out 1 1 -r 2 -c 2
