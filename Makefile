# Makefile for prflab

# Default compiler-optimization level is 0 (no optimizations).
ifndef O
	O = 0
endif

# Optimization is made configurable in case you want to compare your
# optimizations to optimizations performed by the C compiler.

# Compiler and compilation-flags.
CC = gcc
CFLAGS = -march=native -Wall -O$(O) -D __O$(O)__

# That __On__ is a variable for the C preprocessor.
# Compiler will pick different hard-coded constants for baseline
# running times based on optimizaiton levels.

# include libm - library for maths
LIBS = -lm

OBJS = driver.o kernels.o blend.o smooth.o fcyc.o clock.o

all: driver

driver:	$(OBJS) fcyc.h clock.h defs.h config.h
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o driver

clean: 
	-rm -f $(OBJS) driver core *~ *.o
