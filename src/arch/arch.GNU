# ============================================
# GNU Compiler Collection
#
#
# ============================================
LIBS	= 

CFLAGS	= -O3 -ansi -pedantic -s -funroll-loops -fopenmp -DGEO_TUI_
CXXFLAGS= -O3 -ansi -pedantic -s -funroll-loops -fopenmp -DGEO_TUI_
LDFLAGS	= -fopenmp 

CXX	= g++
CC	= gcc
FC	= gfortran

