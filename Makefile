#
# comments...
#

# c++ compiler
#CXX=		g++
#CXX=		icpc

# mpi...
CXX=		mpicxx
CC =		mpicc

CXXFLAGS=	-c -g -Wall -O3 
SOURCES=	main.cpp network_solver.cpp network_info.cpp table.cpp mylapack.cpp 
OBJECTS=	$(SOURCES:.cpp=.o)
LDFLAGS=
LIBS=-L/home/shima86/root_dir/LIB -llapack -lblas #-lgfortran	
EXE=	main


all : $(SOURCES) $(EXE)  

$(EXE) : $(OBJECTS)
	$(CXX) $(LDFLAGS) -o  $(EXE) $(OBJECTS) $(LIBS)

.cpp.o:
	$(CXX)  $(CXXFLAGS) $< -o $@ 

clean : 
	rm -rf *.o *.dat *.txt  $(EXE) *.out.* 
