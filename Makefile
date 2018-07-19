#VER=SEQ
#VER=OpenMP
#VER=MPIOpenMP
VER=OpenACC
#VER=ComplexClass

CXXFLAGS=-std=c++11
LINKFLAGS=-std=c++11

#Sequential version
ifeq ($(VER), SEQ)
    EXE = gppKerSeq.ex
    SRC = gppKerSeq.cpp 
endif

#OpenMP3.5 version
ifeq ($(VER), OpenMP)
    EXE = gppOpenMP3.ex
    SRC = gppOpenMP3.cpp 
    CXXFLAGS += -fopenmp
    LINKFLAGS += -fopenmp
endif

#MPI version
ifeq ($(VER), MPIOpenMP)
    EXE = gppMPIOpenMP.ex
    SRC = gppMPIOpenMP3.cpp 
    CXXFLAGS += -fopenmp
    LINKFLAGS += -fopenmp
endif

#Complex class + gpp version
ifeq ($(VER), OpenACC)
    EXE = gppOpenACC.ex
    SRC = gppOpenACC.cpp 
#    CXXFLAGS += -fopenacc
#    LINKFLAGS += -fopenacc
endif

#Complex class + gpp version
ifeq ($(VER), ComplexClass)
    EXE = gppComplex.ex
    SRC = gppComplex.cpp 
endif

CXX = CC
CXX = CC # scorep-CC
LINK = ${CXX}
CXXFLAGS+= -O3 -g

ifeq ($(VER), OpenACC)
    CXXFLAGS+=-ta=nvidia,cc35,ptxinfo,cuda9.0,pinned -Minfo=accel
    LINKFLAGS+=-ta=nvidia,cc35,ptxinfo,cuda9.0,pinned -Minfo=accel
    #CXXFLAGS+= -h pragma=acc
endif

OBJ = $(SRC:.cpp=.o)

$(EXE): $(OBJ)  
	$(CXX) $(OBJ) -o $(EXE) $(LINKFLAGS)

$(OBJ1): $(SRC) 
	$(CXX) -c $(SRC) $(CXXFLAGS)

clean: 
	rm -f $(OBJ) $(EXE)

