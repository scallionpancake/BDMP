# Use Intel C compiler
CC = icc

# Compile-time flags
CFLAGS = -fast

default: BDMP_sk.c 
	$(CC) $(CFLAGS) -openmp -o paramag BDMP_sk.c utility.c

clean:
	rm paramag
