#VER=SEQ
#VER=OpenMP
VER=MPIOpenMP
#VER=OpenACC
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
endif

#Complex class + gpp version
ifeq ($(VER), OpenACC)
    EXE = gppOpenACC.ex
    SRC = gppOpenACC.cpp 
endif

#Complex class + gpp version
ifeq ($(VER), ComplexClass)
    EXE = gppComplex.ex
    SRC = gppComplex.cpp 
endif

CXX = CC
LINK = ${CXX}
CXXFLAGS+= -O3 -g

ifeq ($(VER), OpenACC)
    CXXFLAGS+=-h pragma=acc
endif

OBJ = $(SRC:.cpp=.o)

$(EXE): $(OBJ)  
	$(CXX) $(OBJ) -o $(EXE) $(LINKFLAGS)

$(OBJ1): $(SRC) 
	$(CXX) -c $(SRC) $(CXXFLAGS)

clean: 
	rm -f $(OBJ) $(EXE)

