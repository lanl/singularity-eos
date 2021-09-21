# © 2021. Triad National Security, LLC. All rights reserved.  This
# program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by Triad
# National Security, LLC for the U.S.  Department of Energy/National
# Nuclear Security Administration. All rights in the program are
# reserved by Triad National Security, LLC, and the U.S. Department of
# Energy/National Nuclear Security Administration. The Government is
# granted for itself and others acting on its behalf a nonexclusive,
# paid-up, irrevocable worldwide license in this material to reproduce,
# prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.

# on snow:
# module purge
# module load intel/19.0.4 hdf5-serial

CC=h5c++

SPINER=../spiner/ports-of-call/portable_arrays.hpp \
	../spiner/ports-of-call/portability.hpp \
	../spiner/databox.hpp \
	../spiner/interpolation.hpp

HDRS=io_eospac.hpp parser.hpp generate_files.hpp ${SPINER}
OBJ=io_eospac.o generate_files.o parser.o
#
CATCH_INCLUDE=-I../catch2
FMATH_INCLUDE=-I../herumi-fmath -I../fast-math
SPINER_INCLUDE=-I../spiner -I../spiner/ports-of-call
#
HDF_LOC=/home/jonahm/hdf5-1.8.20/hdf5
HDF_INCLUDE=$(HDF_LOC)/include
HDF_LIB_LOC=$(HDF_LOC)/lib
HDF_LIB=-L$(HDF_LIB_LOC)
HDF_LINK=-lhdf5_hl -lhdf5 -Wl,-rpath=$(HDF_LIB_LOC)

# snow
# EOSPAC_INCLUDE=-I/usr/projects/data/eos/eospac/eospac6-latest/include/x86_64-avx2-rhel7_toss3.4-4-linux-gnu/intel/19.0.4/
# EOSPAC_LIB_LOC=/usr/projects/data/eos/eospac/eospac6-latest/lib/x86_64-avx2-rhel7_toss3.4-4-linux-gnu/intel/19.0.4/
# darwin
EOSPAC_INCLUDE=-I/usr/projects/data/eos/eospac/eospac6-latest/include/x86_64-avx2-centos7-linux-gnu/gcc/7.4.0/
EOSPAC_LIB_LOC=/usr/projects/data/eos/eospac/eospac6-latest/lib/x86_64-avx2-centos7-linux-gnu/gcc/7.4.0/
# local installation
# EOSPAC_INCLUDE=-I/home/jonahm/eospac6.4.1beta/Source/include/x86_64-avx2-centos7-linux-gnu/gcc/7.4.0/
# EOSPAC_LIB_LOC=/home/jonahm/eospac6.4.1beta/Source/lib/x86_64-avx2-centos7-linux-gnu/gcc/7.4.0/
EOSPAC_LIB=-L${EOSPAC_LIB_LOC}
EOSPAC_LINK=-leospac6 -Wl,-rpath=${EOSPAC_LIB_LOC}

BASE_FLAGS=-std=c++14 -Wall -g -DEOSPAC_WARN
INCLUDE_FLAGS=${CATCH_INCLUDE} ${EOSPAC_INCLUDE} ${FMATH_INCLUDE}
INCLUDE_FLAGS+=${SPINER_INCLUDE}
INCLUDE_FLAGS+=-I${HDF_INCLUDE}

LIB_FLAGS=$(HDF_LIB) $(EOSPAC_LIB)
LINK_FLAGS=$(HDF_LINK) $(EOSPAC_LINK)
CFLAGS=${BASE_FLAGS}

default: sesame2spiner

all: test example

example: sesame2spiner examples/air.dat examples/steel.dat
	./sesame2spiner -p materials.sp5 examples/air.dat examples/steel.dat

test: test.bin
	./test.bin

sesame2spiner: main.o ${OBJ}
	${CXX} ${CFLAGS} ${LIB_FLAGS} -o $@ $^ ${LINK_FLAGS}

test.bin: test.o ${OBJ}
	${CXX} ${CFLAGS} ${LIB_FLAGS} -o $@ $^ ${LINK_FLAGS}

%.o: %.cpp ${HDRS}
	${CXX} ${BASE_FLAGS} ${INCLUDE_FLAGS} -c -o $@ $<

#.PHONY: test clean default example all

clean:
	${RM} test.bin sesame2spiner
	${RM} test.o main.o ${OBJ}
	${RM} sesame_table.sp5 materials.sp5
