HIPCC = hipcc 
CPPFLAGS = -g -O3
HIP_INCLUDE = -I${ROCM_PATH}/include 
LDFLAGS = -L${ROCM_PATH}/lib -lrocsparse -lrocblas
BUILD_DIR = solve


ifeq (${HIP_PLATFORM}, nvcc)
    CPPFLAGS += -arch=compute_86
endif

GPU_SRC = $(wildcard *.cpp)
GPU_OBJ = $(addprefix ${BUILD_DIR}/,$(subst .cpp,.o, $(GPU_SRC)))
GPU_BIN = ${BUILD_DIR}/gpusolver

.PHONY: all clean run solver

all: ${GPU_BIN}

${GPU_BIN}: ${GPU_OBJ}
	${HIPCC} ${LDFLAGS} -o ${GPU_BIN} ${GPU_OBJ}

${BUILD_DIR}/%.o: %.cpp Makefile
	mkdir -p ${BUILD_DIR}
	${HIPCC} ${HIP_INCLUDE} ${CPPFLAGS} -c -o $@ $<  

run: solver
solver:
	${GPU_BIN} $(nd) $(rpd) $(nnzbs) $(block_size)

clean:
	rm -rf ${BUILD_DIR}
