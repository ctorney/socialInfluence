
# Top-Level CmakeLists.txt
cmake_minimum_required(VERSION 2.8)
PROJECT( INF )

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


ADD_SUBDIRECTORY( visuals )

INCLUDE_DIRECTORIES(${INF_SOURCE_DIR}/visuals)

ADD_EXECUTABLE(socialInfluence main)

TARGET_LINK_LIBRARIES(socialInfluence plotGrid ${VTK_LIBRARIES})
