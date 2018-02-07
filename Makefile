#
UNAME = $(shell uname)
ifeq ($(UNAME),Linux)
  BLAS_ROOT = /opt/intel/compilers_and_libraries/linux/mkl
  BLAS_INC_DIR = $(BLAS_ROOT)/include
  BLAS_LIB_DIR = $(BLAS_ROOT)/lib/intel64
  SBLAS_LIBS = -Wl,--start-group $(BLAS_LIB_DIR)/libmkl_intel_lp64.a $(BLAS_LIB_DIR)/libmkl_sequential.a $(BLAS_LIB_DIR)/libmkl_core.a -Wl,--end-group -lpthread -ldl -lm
  TMATRIX_ROOT = /home/stomo/cuda-workspace/Remote/TileMatrix
  COREBLAS_ROOT = /home/stomo/cuda-workspace/Remote/CoreBlas
  CXX =	g++
endif
ifeq ($(UNAME),Darwin)
  BLAS_ROOT = /opt/intel/compilers_and_libraries/mac/mkl
  BLAS_INC_DIR = $(BLAS_ROOT)/include
  BLAS_LIB_DIR = $(BLAS_ROOT)/lib
  SBLAS_LIBS = $(BLAS_LIB_DIR)/libmkl_intel_lp64.a $(BLAS_LIB_DIR)/libmkl_sequential.a $(BLAS_LIB_DIR)/libmkl_core.a -lpthread -ldl -lm
  TMATRIX_ROOT = /Users/stomo/WorkSpace/TileAlgorithm/TileMatrix
  COREBLAS_ROOT = /Users/stomo/WorkSpace/TileAlgorithm/CoreBlas
  CXX =	/usr/local/bin/g++ 
endif

#
PLASMA_ROOT = /opt/plasma-17.1
PLASMA_INC_DIR = $(PLASMA_ROOT)/include
PLASMA_LIB_DIR = $(PLASMA_ROOT)/lib
PLASMA_LIBS = -lcoreblas -lplasma
#
TMATRIX_INC_DIR = $(TMATRIX_ROOT)
TMATRIX_LIB_DIR = $(TMATRIX_ROOT)
TMATRIX_LIBS = -lTileMatrix
#
COREBLAS_INC_DIR = $(COREBLAS_ROOT)
COREBLAS_LIB_DIR = $(COREBLAS_ROOT)
COREBLAS_LIBS = -lCoreBlasTile
#
CXXFLAGS =	-fopenmp -m64 -O2 -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR) -I$(TMATRIX_INC_DIR) -I$(COREBLAS_INC_DIR)
#
TDOBJS =		TileTRD.o TRD.o
TTOBJS =		TileTRD.o TRD_Task.o
SDOBJS =		TileTRD.o TRD_Sym.o
STOBJS =		TileTRD.o TRD_Sym_Task.o

all: ST SD TT TD

ST:	$(STOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(STOBJS) \
				-L$(TMATRIX_LIB_DIR) $(TMATRIX_LIBS) \
				-L$(COREBLAS_LIB_DIR) $(COREBLAS_LIBS) \
				-L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) \
				-L$(BLAS_LIB_DIR) $(SBLAS_LIBS)

SD:	$(SDOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(SDOBJS) \
				-L$(TMATRIX_LIB_DIR) $(TMATRIX_LIBS) \
				-L$(COREBLAS_LIB_DIR) $(COREBLAS_LIBS) \
				-L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) \
				-L$(BLAS_LIB_DIR) $(SBLAS_LIBS)

TT:	$(TTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(TTOBJS) \
				-L$(TMATRIX_LIB_DIR) $(TMATRIX_LIBS) \
				-L$(COREBLAS_LIB_DIR) $(COREBLAS_LIBS) \
				-L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) \
				-L$(BLAS_LIB_DIR) $(SBLAS_LIBS)

TD:	$(TDOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(TDOBJS) \
				-L$(TMATRIX_LIB_DIR) $(TMATRIX_LIBS) \
				-L$(COREBLAS_LIB_DIR) $(COREBLAS_LIBS) \
				-L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) \
				-L$(BLAS_LIB_DIR) $(SBLAS_LIBS)

clean:
	rm -f $(TDOBJS) $(TTOBJS) $(SDOBJS) $(STOBJS) ST SD TT TD
