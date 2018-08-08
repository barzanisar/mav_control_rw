cmake_minimum_required(VERSION 2.8.3)
project(mav_modified_rodrigues_nmpc)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")

# On os x we need to use openblas for the lapacke library
find_package(OpenBLAS REQUIRED)
if (Openblas_FOUND)
  message("-- Found Openblas library: ${Openblas_LIBRARIES}")
  message("-- Found Openblas include dirs: ${Openblas_INCLUDE_DIRS}")
  include_directories(${Openblas_INCLUDE_DIRS})
  set(lapacke_LIBRARIES ${Openblas_LIBRARIES})
else (Openblas_FOUND)
  message("-- Openblas not found, trying to use lapacke")
  find_package(LAPACK REQUIRED)
  set(lapacke_LIBRARIES lapacke)
endif (Openblas_FOUND)


find_package(catkin_simple REQUIRED)
catkin_simple(ALL_DEPS_REQUIRED)

if(APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

find_package(Eigen REQUIRED)
include_directories(${Eigen_INCLUDE_DIRS})

#############
# LIBRARIES #
#############
cs_add_library(mav_modified_rodrigues_nmpc_lib
  src/modified_rodrigues_nmpc.cc
  solver/qpoases/SRC/Bounds.cpp
  solver/qpoases/SRC/Constraints.cpp
  solver/qpoases/SRC/CyclingManager.cpp
  solver/qpoases/SRC/Indexlist.cpp
  solver/qpoases/SRC/MessageHandling.cpp
  solver/qpoases/SRC/QProblem.cpp
  solver/qpoases/SRC/QProblemB.cpp
  solver/qpoases/SRC/SubjectTo.cpp
  solver/qpoases/SRC/Utils.cpp
  solver/qpoases/SRC/EXTRAS/SolutionAnalysis.cpp
  solver/acado_qpoases_interface.cpp
  solver/acado_integrator.c
  solver/acado_solver.c
  solver/acado_auxiliary_functions.c
  )

target_include_directories(mav_modified_rodrigues_nmpc_lib PUBLIC
	solver/
	solver/qpoases
	solver/qpoases/INCLUDE
	solver/qpoases/SRC
	)

############
# BINARIES #
############
cs_add_executable(modified_rodrigues_nmpc_node
  src/modified_rodrigues_nmpc_node.cc
)



target_link_libraries(mav_modified_rodrigues_nmpc_lib ${lapacke_LIBRARIES})

target_link_libraries(modified_rodrigues_nmpc_node mav_modified_rodrigues_nmpc_lib)

add_dependencies(modified_rodrigues_nmpc_node ${${PROJECT_NAME}_EXPORTED_TARGETS})

cs_add_executable(ref_trajectory_node
  src/ref_trajectory_node.cpp
)

#target_link_libraries(ref_trajectory_node mav_modified_rodrigues_nmpc_lib)

##########
# EXPORT #
##########
cs_install()
cs_export()
