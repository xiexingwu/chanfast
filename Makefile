SDIR = ./src
IDIR = ./include
ODIR = ./src/obj
LOGS = ./logs

CC = mpicc
CFLAGS = -I$(IDIR)
LIBS = -lm -lfftw3
WARNS = -Wno-discarded-qualifiers -Wno-int-to-pointer-cast -Wno-incompatible-pointer-types

_DEP = check.h logging.h parser.h globals_sim.h fourth.h rk3.h
_OBJ = main.o logging.o mpi_plan.o tensor.o

TARGET = a.out

DEP = $(patsubst %,$(IDIR)/%,$(_DEP))
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(DEFINES) $(INCLUDES)

$(ODIR)/%.o: $(SDIR)/%.c $(DEP)
	$(CC) -c -o $@ $< $(CFLAGS) $(DEFINES) $(INCLUDES) $(WARNS)

.PHONY: clean clobber all
clean:
	rm -f $(ODIR)/*.o $(LOGS)/*

clobber: clean
	rm -f $(TARGET)

all: clobber $(TARGET)