#------------------------------------------------------------------------------
# © 2021 - 2022. Triad National Security, LLC. All rights reserved.  This
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
#------------------------------------------------------------------------------


#USE_HDF5=yes
USE_HDF5=no
USE_KOKKOS=no
USE_EOSPAC=no
INVERT_AT_SETUP=no

CXXFLAGS   := -O3 -DNDEBUG
CXX        := g++
# CXX        := mpic++
#CXXFLAGS   := -g -DEOSPAC_WARN
#CXXFLAGS   := -g -DEOSPAC_WARN -fopenmp
LDFLAGS    := #-lstdc++ -lm
LIBFLAGS   :=
BUILD      := ./build
OBJ_DIR    := $(BUILD)/objects
APP_DIR    := $(BUILD)/apps
DATA_DIR   := ../data
TARGET_PTE := test_pte
TARGET_SPN := test_eos_spiner
TARGET_CMP := compare_to_eospac
TARGET_PROF := profile_eos
TARGET_DATA := compare_on_data
TARGET_DATASET := compare_on_dataset
INCLUDE    := -I../utils/eigen/ -I../utils/catch2/ -I../utils/nlohmann -I../utils/ports-of-call -I../utils/spiner -I../utils/herumi-fmath -I../utils/fast-math

HDF_LOC=/home/jonahm/hdf5-1.8.20/hdf5
#HDF_INCLUDE=$(HDF_LOC)/include
HDF_INCLUDE=/usr/include/hdf5/openmpi/
#HDF_LIB_LOC=$(HDF_LOC)/lib
HDF_LIB_LOC= /usr/lib/x86_64-linux-gnu/hdf5/openmpi/
HDF_LIB=-L$(HDF_LIB_LOC)
HDF_LINK=-lhdf5_hl -lhdf5 -Wl,-rpath=$(HDF_LIB_LOC)

# for snow
# EOSPAC_INCLUDE=-I/usr/projects/data/eos/eospac/eospac6-latest/include/x86_64-avx2-rhel7_toss3.4-4-linux-gnu/intel/19.0.4/
# EOSPAC_LIB_LOC=/usr/projects/data/eos/eospac/eospac6-latest/lib/x86_64-avx2-rhel7_toss3.4-4-linux-gnu/intel/19.0.4/
# for darwin
EOSPAC_INCLUDE=-I/usr/projects/data/eos/eospac/eospac6-latest/include/x86_64-avx2-centos7-linux-gnu/gcc/7.4.0
EOSPAC_LIB_LOC=/usr/projects/data/eos/eospac/eospac6-latest/lib/x86_64-avx2-centos7-linux-gnu/gcc/7.4.0
EOSPAC_LIB=-L${EOSPAC_LIB_LOC}
EOSPAC_LINK=-leospac6 -Wl,-rpath=${EOSPAC_LIB_LOC}

KOKKOS_PATH = /home/jonahm/kokkos
KOKKOS_CXX_STANDARD=c++14
KOKKOS_DEVICES = "Cuda"
KOKKOS_CUDA_OPTIONS = "enable_lambda;rdc"
KOKKOS_ARCH = "Volta70"

ifeq ($(USE_KOKKOS), yes)
	SRC := $(wildcard ../eos/*.cpp)
else
	SRC := $(wildcard ../eos/*.cpp) $(wildcard ../mix/*.cpp)
endif

ifeq ($(USE_HDF5),yes)
	INCLUDE += -I$(HDF_INCLUDE)
	LIBFLAGS += $(HDF_LIB)
	LDFLAGS += $(HDF_LINK)
	CXXFLAGS := $(CXXFLAGS) -DSPINER_USE_HDF
else
	CXXFLAGS := $(CXXFLAGS)
endif

ifeq ($(USE_EOSPAC),yes)
	SRC += ../utils/sesame2spiner/io_eospac.cpp
	CXXFLAGS := $(CXXFLAGS) -DUSE_EOSPAC -DEOSPAC_WARN
	INCLUDE += $(EOSPAC_INCLUDE)
	LIBFLAGS += $(EOSPAC_LIB)
	LDFLAGS += $(EOSPAC_LINK)
	ifeq ($(INVERT_AT_SETUP),yes)
		CXXFLAGS := $(CXXFLAGS) -DTEST_EOSPAC_INVERT_AT_SETUP
	endif
endif

# kokkos stuff
ifeq ($(USE_KOKKOS),yes)
	CXXFLAGS += -DPORTABILITY_STRATEGY_KOKKOS --expt-relaxed-constexpr
	CXX := $(KOKKOS_PATH)/bin/nvcc_wrapper
	include $(KOKKOS_PATH)/Makefile.kokkos
else
	KOKKOS_CXXFLAGS=
	KOKKOS_CPPFLAGS=
	KOKKOS_LINK_DEPENDS=
	KOKKOS_CPP_DEPENDS=
	KOKKOS_LDFLAGS=
	CXXFLAGS += -std=c++14
endif

TMP_SRC := $(SRC)
MAINS := $(wildcard ../utils/*/main.cpp) \
	$(wildcard ../utils/*/test.cpp)
SRC := $(filter-out $(MAINS), $(TMP_SRC))
OBJECTS := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

default: test

test: test_spiner test_pte

all: pte spiner compare profile data dataset

test_pte: pte
	$(APP_DIR)/$(TARGET_PTE)

test_spiner: spiner
	$(APP_DIR)/$(TARGET_SPN)

profile: profile_eos plot_timings.py
	$(APP_DIR)/$(TARGET_PROF) 100 8 16 32 64 128 256
	python plot_timings.py IdealGas_timing.dat -s ../figures/ideal-gas-timings

pte: build $(APP_DIR)/$(TARGET_PTE)

spiner: build $(APP_DIR)/$(TARGET_SPN)

compare: build $(APP_DIR)/$(TARGET_CMP)

profile_eos: build $(APP_DIR)/$(TARGET_PROF)

data: build $(APP_DIR)/$(TARGET_DATA)

dataset: build $(APP_DIR)/$(TARGET_DATASET)

$(OBJ_DIR)/%.o: %.cpp $(KOKKOS_CPP_DEPENDS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) -c $< -o $@ 

$(APP_DIR)/$(TARGET_PTE): $(OBJECTS) $(OBJ_DIR)/$(TARGET_PTE).o
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBFLAGS) -o $@ $^ $(LDFLAGS)

$(APP_DIR)/$(TARGET_SPN): $(OBJECTS) $(OBJ_DIR)/$(TARGET_SPN).o $(KOKKOS_LINK_DEPENDS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBFLAGS) -o $@ $^ $(KOKKOS_LDFLAGS) $(LDFLAGS)

$(APP_DIR)/$(TARGET_CMP): $(OBJECTS) $(OBJ_DIR)/$(TARGET_CMP).o $(KOKKOS_LINK_DEPENDS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBFLAGS) -o $@ $^ $(KOKKOS_LDFLAGS) $(LDFLAGS)

$(APP_DIR)/$(TARGET_PROF): $(OBJECTS) $(OBJ_DIR)/$(TARGET_PROF).o $(KOKKOS_LINK_DEPENDS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBFLAGS) -o $@ $^ $(KOKKOS_LDFLAGS) $(LDFLAGS)

$(APP_DIR)/$(TARGET_DATA): $(OBJECTS) $(OBJ_DIR)/$(TARGET_DATA).o $(KOKKOS_LINK_DEPENDS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBFLAGS) -o $@ $^ $(KOKKOS_LDFLAGS) $(LDFLAGS)

$(APP_DIR)/$(TARGET_DATASET): $(OBJECTS) $(OBJ_DIR)/$(TARGET_DATASET).o $(KOKKOS_LINK_DEPENDS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBFLAGS) -o $@ $^ $(KOKKOS_LDFLAGS) $(LDFLAGS)

.PHONY: all build clean debug release pte spiner test default test_spiner test_pte profile_eos profile compare data

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O2
release: all

clean:
	@rm -rvf $(BUILD)
	@rm -rvf *.o *.a
	@rm -rvf KokkosCore_config.*
