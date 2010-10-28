
include mak.inc

LIBS	= 

CXXINC	=
INC	=
DEFS	= 
CFLAGS	= -g -Wno-deprecated
LDFLAGS	= 

CXX	= g++
CC	= gcc
RM	= rm -f

.C.o: ; $(CXX) -c $(CFLAGS) $*.C
.c.o: ; $(CC) -c $(CFLAGS) $*.c

all: $(EXEC)

$(EXEC): $(OBJS) 
	$(CXX) $(LDFLAGS) -o $@ $(INC) $(OBJS) $(LIBS)

.PHONY : clean
clean:
	-$(RM) *.o *~ $(EXEC)

depend:
	makedepend $(INCLUDES) *.C

