SDIR = ./src
IDIR = ./include
ODIR = ./src/obj

CC = mpicc
CFLAGS = -I$(IDIR)
LIBS = -lm

_DEP = check.h parser.h globals_sim.h 
_OBJ = mpi_plan.o main.o 

TARGET = a.out

DEP = $(patsubst %,$(IDIR)/%,$(_DEP))
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(DEFINES)

$(ODIR)/%.o: $(SDIR)/%.c $(DEP)
	$(CC) -c -o $@ $< $(CFLAGS) $(DEFINES)

.PHONY: clean clobber all
clean:
	rm -f $(ODIR)/*.o 

clobber: clean
	rm -f $(TARGET)

all: clobber $(TARGET)