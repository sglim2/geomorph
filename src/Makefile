
include arch/mak.inc
include arch/arch.$(ARCH)

OBJS	= $(OBJS1) $(OBJS2)
RM	= rm -f

.cpp.o: ; $(CXX) -c $(CXXFLAGS) $*.cpp
.c.o: ; $(CC) -c $(CFLAGS) $*.c
.f.o: ; $(FC) -c $(FFLAGS) $*.f

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $(INC) $(OBJS) $(LIBS) 

.PHONY : clean
clean:
	-$(RM) *.o *~ $(EXEC)

arch. : 
# Complain when ARCH is missing 
	@echo " "
	@echo " You need to compile with an ARCH file"
	@echo " make ARCH=[]"
	@echo " Options include:"
	@echo " [ GNU ]"
	@echo " [ OSX ]"
	@echo " [ INTEL]"
	@echo " [ PGI ]"
	@echo " [ DEBUG ]"
	@/dev/null/abort >& /dev/null

