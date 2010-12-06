# ============================================
# PGI Compiler Suite; no GPU support.
#
#
# ============================================
LIBS	= 

CFLAGS	= -O3 -ansi -pedantic -s -funroll-loops -fopenmp
CXXFLAGS= -O3 -ansi -pedantic -s -funroll-loops -fopenmp
LDFLAGS	= -fopenmp

CXX	= g++
CC	= gcc
FC	= gfortran

