# Compiler settings
CC = gcc
FC = gfortran
CFLAGS = -c -g
FFLAGS = -g                  # Fortran 也要加上 -g
LDFLAGS = -lX11
INCLUDES = -I.        # current directory for headers

# File names
C_SRC = pslib.c tinyspline.c parson.c
C_OBJ = pslib.o tinyspline.o parson.o
F_SRC = spline_demo.f
EXE = spline_demo

# Default target
all: $(EXE)

# Compile C sources to object files
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@

# Link Fortran and C objects to create the executable
$(EXE): $(F_SRC) $(C_OBJ)
	$(FC) $(FFLAGS) $(F_SRC) $(C_OBJ) -o $(EXE) $(LDFLAGS)

# Clean up
clean:
	rm -f $(C_OBJ) $(EXE)
