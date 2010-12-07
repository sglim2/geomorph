# ============================================
# PGI Compiler Suite; no GPU support.
#
#
# ============================================
LIBS	= 

OBJS2	= dummy.o
CFLAGS	= -O3 -ansi -pedantic -s -funroll-loops -fopenmp
CXXFLAGS= -O3 -ansi -pedantic -s -funroll-loops -fopenmp
LDFLAGS	= -fopenmp

CXX	= g++
CC	= gcc
FC	= gfortran

