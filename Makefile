# Modern OP2 build for BookLeaf.
#
# Source files remain under src/.  The OP2 translator writes every backend
# product under generated/bookleaf/, so no generated source is checked in.

include $(OP2_INSTALL_PATH)/../makefiles/common.mk

APP_NAME := bookleaf

APP_SRC_PRE_KERNEL := \
  src/utils/op2_bookleaf_api.F90 \
  src/utils/data.F90 \
  src/utils/bookleaf_consts.F90 \
  src/common_kernels.F90 \
  src/utils/error.F90 \
  src/utils/timers.F90 \
  src/utils/geometry_kernels.F90 \
  src/utils/init_kernels.F90 \
  src/eos/eos.F90 \
  src/eos/getpc_kernels.F90 \
  src/hydro/gethg_kernels.F90 \
  src/hydro/gethg.F90 \
  src/hydro/getsp_kernels.F90 \
  src/hydro/getsp.F90 \
  src/hydro/getq_kernels.F90 \
  src/getacc_kernels.F90 \
  src/getdt_kernels.F90 \
  src/getein_kernels.F90 \
  src/getforce_kernels.F90 \
  src/lagstep_kernels.F90 \
  src/ale/ale_advectors_kernels.F90 \
  src/ale/ale_advect_kernels.F90 \
  src/ale/ale_getfvol_kernels.F90 \
  src/ale/ale_getmesh_kernels.F90 \
  src/ale/ale_update_kernels.F90 \
  src/io/write_kernels.F90

APP_SRC_POST_KERNEL := \
  src/utils/utilities.F90 \
  src/utils/geometry.F90 \
  src/eos/getpc.F90 \
  src/hydro/getq.F90 \
  src/getacc.F90 \
  src/getdt.F90 \
  src/getein.F90 \
  src/getforce.F90 \
  src/lagstep.F90 \
  src/ale/ale_advectors.F90 \
  src/ale/ale_advect.F90 \
  src/ale/ale_getfvol.F90 \
  src/ale/ale_getmesh.F90 \
  src/ale/ale_update.F90 \
  src/ale/alestep.F90 \
  src/utils/init.F90 \
  src/utils/mesh.F90 \
  src/io/banner.F90 \
  src/io/read.F90 \
  src/io/write.F90 \
  src/utils/halt.F90 \
  src/hydro.F90 \
  src/main.F90

APP_EXTRA_TRANSLATOR_FLAGS := -D OP2_TRANSLATOR -t seq -t openmp -t c_cuda -t c_hip --consts-module src/utils/bookleaf_consts.F90

OP2_LIBS_WITH_HDF5 := true

# translator-v2 writes all transformed program sources into one directory,
# whereas BookLeaf keeps its source modules in subdirectories.  Its module
# ordering and flattened generated-source paths are therefore application
# concerns; do not require a BookLeaf exception in OP2-Common's f_app.mk.
TRANSLATOR ?= $(ROOT_DIR)/translator-v2/op2-translator.sh -v
BOOKLEAF_SRC := $(APP_SRC_PRE_KERNEL) $(APP_SRC_POST_KERNEL)
BOOKLEAF_GEN_DIR := generated/$(APP_NAME)
BOOKLEAF_GEN_STAMP := $(BOOKLEAF_GEN_DIR)/.generated
BOOKLEAF_GEN_PRE := $(addprefix $(BOOKLEAF_GEN_DIR)/,$(notdir $(APP_SRC_PRE_KERNEL)))
BOOKLEAF_GEN_POST := $(addprefix $(BOOKLEAF_GEN_DIR)/,$(notdir $(APP_SRC_POST_KERNEL)))
BOOKLEAF_CPU_FLAGS := $(OMP_FFLAGS)

BOOKLEAF_MESHGEN_BUILD_DIR := build-meshgen
BOOKLEAF_MESHGEN_MOD_DIR := $(BOOKLEAF_MESHGEN_BUILD_DIR)/mod
BOOKLEAF_MESHGEN_SRC := \
  src/utils/op2_bookleaf_api.F90 \
  src/utils/data_gen.F90 \
  src/utils/bookleaf_consts.F90 \
  src/common_kernels.F90 \
  src/utils/error.F90 \
  src/utils/timers.F90 \
  src/utils/geometry_kernels.F90 \
  src/utils/init_kernels.F90 \
  src/utils/utilities.F90 \
  src/utils/geometry.F90 \
  src/utils/mesh.F90 \
  src/io/write_kernels.F90 \
  src/io/write.F90 \
  src/utils/halt.F90 \
  src/io/read.F90 \
  src/utils/init_gen.F90 \
  src/io/banner.F90 \
  src/main_gen.F90

include $(OP2_INSTALL_PATH)/../makefiles/lib_helpers.mk

.PHONY: all clean generate \
  bookleaf_seq bookleaf_genseq bookleaf_openmp bookleaf_c_cuda bookleaf_c_hip \
  bookleaf_mpi_seq bookleaf_mpi_genseq bookleaf_mpi_openmp \
  bookleaf_mpi_c_cuda bookleaf_mpi_c_hip bookleaf_meshgen
.PRECIOUS: $(BOOKLEAF_GEN_STAMP) mod/%

BOOKLEAF_BASE_VARIANTS := seq genseq openmp c_cuda c_hip
BOOKLEAF_ALL_VARIANTS := $(BOOKLEAF_BASE_VARIANTS) \
  $(addprefix mpi_,$(BOOKLEAF_BASE_VARIANTS))
BOOKLEAF_BUILDABLE_VARIANTS :=

ifeq ($(HAVE_F),true)
  BOOKLEAF_BUILDABLE_VARIANTS += seq genseq
  ifeq ($(F_HAS_OMP),true)
    BOOKLEAF_BUILDABLE_VARIANTS += openmp
  endif
  ifeq ($(HAVE_C_CUDA),true)
    BOOKLEAF_BUILDABLE_VARIANTS += c_cuda
  endif
  ifeq ($(HAVE_C_HIP),true)
    BOOKLEAF_BUILDABLE_VARIANTS += c_hip
  endif
  ifeq ($(HAVE_MPI_F),true)
    BOOKLEAF_BUILDABLE_VARIANTS += $(addprefix mpi_,$(BOOKLEAF_BUILDABLE_VARIANTS))
  endif
endif

ifeq ($(OP2_LIBS_WITH_HDF5),true)
  ifneq ($(HAVE_HDF5_SEQ),true)
    BOOKLEAF_BUILDABLE_VARIANTS := $(filter mpi_%,$(BOOKLEAF_BUILDABLE_VARIANTS))
  endif
  ifneq ($(HAVE_HDF5_PAR),true)
    BOOKLEAF_BUILDABLE_VARIANTS := $(filter-out mpi_%,$(BOOKLEAF_BUILDABLE_VARIANTS))
  endif
endif

all: $(addprefix bookleaf_,$(BOOKLEAF_BUILDABLE_VARIANTS))

$(BOOKLEAF_GEN_STAMP): $(BOOKLEAF_SRC)
	@mkdir -p $(BOOKLEAF_GEN_DIR)
	$(TRANSLATOR) $(APP_EXTRA_FLAGS) $(APP_EXTRA_TRANSLATOR_FLAGS) $^ -o $(BOOKLEAF_GEN_DIR)
	@touch $@

generate: $(BOOKLEAF_GEN_STAMP)

mod/%:
	@mkdir -p $@

$(BOOKLEAF_MESHGEN_MOD_DIR):
	@mkdir -p $@

bookleaf_meshgen: $(BOOKLEAF_MESHGEN_SRC) | $(BOOKLEAF_MESHGEN_MOD_DIR)
	$(FC) $(FFLAGS) $(BOOKLEAF_CPU_FLAGS) $(APP_EXTRA_FLAGS) $(F_MOD_OUT_OPT)$(BOOKLEAF_MESHGEN_MOD_DIR) $(OP2_MOD) \
		$(BOOKLEAF_MESHGEN_SRC) $(OP2_LIB_FOR_SEQ) $(CXXLINK) -o $@

# Developer build: no translator invocation and no generated sources.
bookleaf_seq: $(BOOKLEAF_SRC) | mod/$(APP_NAME)/seq
	$(FC) $(FFLAGS) $(BOOKLEAF_CPU_FLAGS) $(APP_EXTRA_FLAGS) $(F_MOD_OUT_OPT)$| $(OP2_MOD) \
		$(BOOKLEAF_SRC) $(OP2_LIB_FOR_SEQ) $(CXXLINK) -o $@

bookleaf_genseq: $(BOOKLEAF_GEN_STAMP) | mod/$(APP_NAME)/genseq
	$(FC) $(FFLAGS) $(BOOKLEAF_CPU_FLAGS) $(APP_EXTRA_FLAGS) $(F_MOD_OUT_OPT)$| $(OP2_MOD) \
		$(BOOKLEAF_GEN_DIR)/seq/op2_consts.F90 $(BOOKLEAF_GEN_PRE) \
		$(BOOKLEAF_GEN_DIR)/seq/op2_kernels.F90 \
		$(BOOKLEAF_GEN_POST) $(OP2_LIB_FOR_SEQ) $(CXXLINK) -o $@

bookleaf_openmp: $(BOOKLEAF_GEN_STAMP) | mod/$(APP_NAME)/openmp
	$(FC) $(FFLAGS) $(BOOKLEAF_CPU_FLAGS) $(APP_EXTRA_FLAGS) $(F_MOD_OUT_OPT)$| $(OP2_MOD) \
		$(BOOKLEAF_GEN_DIR)/openmp/op2_consts.F90 $(BOOKLEAF_GEN_PRE) \
		$(BOOKLEAF_GEN_DIR)/openmp/op2_kernels.F90 \
		$(BOOKLEAF_GEN_POST) $(OP2_LIB_FOR_OPENMP) $(CXXLINK) -o $@

bookleaf_mpi_seq: $(BOOKLEAF_SRC) | mod/$(APP_NAME)/mpi_seq
	$(MPIFC) $(FFLAGS) $(BOOKLEAF_CPU_FLAGS) $(APP_EXTRA_FLAGS) $(F_MOD_OUT_OPT)$| $(OP2_MOD) \
		$(BOOKLEAF_SRC) $(OP2_LIB_FOR_MPI) $(CXXLINK) -o $@

bookleaf_mpi_genseq: $(BOOKLEAF_GEN_STAMP) | mod/$(APP_NAME)/mpi_genseq
	$(MPIFC) $(FFLAGS) $(BOOKLEAF_CPU_FLAGS) $(APP_EXTRA_FLAGS) $(F_MOD_OUT_OPT)$| $(OP2_MOD) \
		$(BOOKLEAF_GEN_DIR)/seq/op2_consts.F90 $(BOOKLEAF_GEN_PRE) \
		$(BOOKLEAF_GEN_DIR)/seq/op2_kernels.F90 \
		$(BOOKLEAF_GEN_POST) $(OP2_LIB_FOR_MPI) $(CXXLINK) -o $@

bookleaf_mpi_openmp: $(BOOKLEAF_GEN_STAMP) | mod/$(APP_NAME)/mpi_openmp
	$(MPIFC) $(FFLAGS) $(BOOKLEAF_CPU_FLAGS) $(APP_EXTRA_FLAGS) $(F_MOD_OUT_OPT)$| $(OP2_MOD) \
		$(BOOKLEAF_GEN_DIR)/openmp/op2_consts.F90 $(BOOKLEAF_GEN_PRE) \
		$(BOOKLEAF_GEN_DIR)/openmp/op2_kernels.F90 \
		$(BOOKLEAF_GEN_POST) $(OP2_LIB_FOR_MPI) $(CXXLINK) -o $@

$(BOOKLEAF_GEN_DIR)/c_cuda/op2_kernels.o: $(BOOKLEAF_GEN_STAMP)
	$(NVCC) $(NVCCFLAGS) $(OP2_INC) $(APP_EXTRA_FLAGS) -DOP2_CUDA \
		-c $(BOOKLEAF_GEN_DIR)/c_cuda/op2_kernels_aux1.cu -o $@

bookleaf_c_cuda: $(BOOKLEAF_GEN_STAMP) $(BOOKLEAF_GEN_DIR)/c_cuda/op2_kernels.o | mod/$(APP_NAME)/c_cuda
	$(FC) $(FFLAGS) $(BOOKLEAF_CPU_FLAGS) $(CUDA_FFLAGS) $(APP_EXTRA_FLAGS) $(F_MOD_OUT_OPT)$| $(OP2_MOD) \
		$(BOOKLEAF_GEN_DIR)/c_cuda/op2_consts.F90 $(BOOKLEAF_GEN_PRE) \
		$(BOOKLEAF_GEN_DIR)/c_cuda/op2_kernels.F90 \
		$(BOOKLEAF_GEN_POST) $(BOOKLEAF_GEN_DIR)/c_cuda/op2_kernels.o \
		$(OP2_LIB_FOR_CUDA) $(CXXLINK) $(CUDA_LIB) -o $@

bookleaf_mpi_c_cuda: $(BOOKLEAF_GEN_STAMP) $(BOOKLEAF_GEN_DIR)/c_cuda/op2_kernels.o | mod/$(APP_NAME)/mpi_c_cuda
	$(MPIFC) $(FFLAGS) $(BOOKLEAF_CPU_FLAGS) $(CUDA_FFLAGS) $(APP_EXTRA_FLAGS) $(F_MOD_OUT_OPT)$| $(OP2_MOD) \
		$(BOOKLEAF_GEN_DIR)/c_cuda/op2_consts.F90 $(BOOKLEAF_GEN_PRE) \
		$(BOOKLEAF_GEN_DIR)/c_cuda/op2_kernels.F90 \
		$(BOOKLEAF_GEN_POST) $(BOOKLEAF_GEN_DIR)/c_cuda/op2_kernels.o \
		$(OP2_LIB_FOR_MPI_CUDA) $(CXXLINK) $(CUDA_LIB) -o $@

$(BOOKLEAF_GEN_DIR)/c_hip/op2_kernels.o: $(BOOKLEAF_GEN_STAMP)
	$(HIPCC) $(HIPCCFLAGS) $(OP2_INC) $(APP_EXTRA_FLAGS) -DOP2_HIP \
		-c $(BOOKLEAF_GEN_DIR)/c_hip/op2_kernels_aux1.cu -o $@

bookleaf_c_hip: $(BOOKLEAF_GEN_STAMP) $(BOOKLEAF_GEN_DIR)/c_hip/op2_kernels.o | mod/$(APP_NAME)/c_hip
	$(FC) $(FFLAGS) $(BOOKLEAF_CPU_FLAGS) $(HIP_FFLAGS) $(APP_EXTRA_FLAGS) $(F_MOD_OUT_OPT)$| $(OP2_MOD) \
		$(BOOKLEAF_GEN_DIR)/c_hip/op2_consts.F90 $(BOOKLEAF_GEN_PRE) \
		$(BOOKLEAF_GEN_DIR)/c_hip/op2_kernels.F90 \
		$(BOOKLEAF_GEN_POST) $(BOOKLEAF_GEN_DIR)/c_hip/op2_kernels.o \
		$(OP2_LIB_FOR_HIP) $(CXXLINK) $(HIP_LIB) -o $@

bookleaf_mpi_c_hip: $(BOOKLEAF_GEN_STAMP) $(BOOKLEAF_GEN_DIR)/c_hip/op2_kernels.o | mod/$(APP_NAME)/mpi_c_hip
	$(MPIFC) $(FFLAGS) $(BOOKLEAF_CPU_FLAGS) $(HIP_FFLAGS) $(APP_EXTRA_FLAGS) $(F_MOD_OUT_OPT)$| $(OP2_MOD) \
		$(BOOKLEAF_GEN_DIR)/c_hip/op2_consts.F90 $(BOOKLEAF_GEN_PRE) \
		$(BOOKLEAF_GEN_DIR)/c_hip/op2_kernels.F90 \
		$(BOOKLEAF_GEN_POST) $(BOOKLEAF_GEN_DIR)/c_hip/op2_kernels.o \
		$(OP2_LIB_FOR_MPI_HIP) $(CXXLINK) $(HIP_LIB) -o $@

clean:
	-$(RM) $(addprefix bookleaf_,$(BOOKLEAF_ALL_VARIANTS)) bookleaf_meshgen
	-$(RM) -rf generated mod $(BOOKLEAF_MESHGEN_BUILD_DIR)
