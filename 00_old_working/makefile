# Use Intel C compiler
CC = icc

# Compile-time flags
CFLAGS = -fast

default: buckling_2D_srfrpl_openMP.c 
	$(CC) $(CFLAGS) -openmp -o b2 buckling_2D_srfrpl_openMP.c utility.c

clean:
	rm b2
