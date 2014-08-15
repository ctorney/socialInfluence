
# Top-Level CmakeLists.txt
cmake_minimum_required(VERSION 2.8)
PROJECT( INF )

find_package(VTK REQUIRED)
find_package(CUDA REQUIRED)
include(${VTK_USE_FILE})

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-O3 )

ADD_SUBDIRECTORY( visuals )

INCLUDE_DIRECTORIES(${INF_SOURCE_DIR}/include)

CUDA_ADD_EXECUTABLE(socialInfluence main switcherKernel.cu cnpy.cpp)

TARGET_LINK_LIBRARIES(socialInfluence plotGrid ${VTK_LIBRARIES})
