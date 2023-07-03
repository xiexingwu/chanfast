cmake_minimum_required(VERSION 3.14)
project(chanfast CXX)

# MPI
find_package(MPI REQUIRED)
message(STATUS "mpi found? ${MPI_CXX_FOUND}. -I at ${MPI_CXX_INCLUDE_DIRS}. -L at ${MPI_CXX_LIBRARIES}")

# FFTW
list(APPEND CMAKE_MODULE_PATH 
  ${CMAKE_CURRENT_LIST_DIR}/cmake/findfftw
)
find_package(FFTW REQUIRED)
message(STATUS "fftw found? ${FFTW_FOUND}. -I at ${FFTW_INCLUDE_DIRS}. -L at ${FFTW_LIBRARIES}")

include(FetchContent)
# spdlog
FetchContent_Declare(
  spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG        v1.11.0
)
FetchContent_MakeAvailable(spdlog)

# chanfast
set(chanfast_INCLUDES 
  ${chanfast_SOURCE_DIR}/include
  ${chanfast_SOURCE_DIR}/include/chanfast/consts
  ${chanfast_SOURCE_DIR}/include/chanfast/params
  ${MPI_CXX_INCLUDE_DIRS}
  ${spdlog_SOURCE_DIR}/include
)
message(STATUS "includes ${chanfast_INCLUDES}")

set(chanfast_LIBS
  ${MPI_CXX_LIBRARIES}
  spdlog
)

# src
add_subdirectory(src)

# tests
add_subdirectory(tests)

# main
add_executable(chanfast main.cpp)
target_include_directories(chanfast PRIVATE ${chanfast_INCLUDES})
target_link_libraries(chanfast PRIVATE ${chanfast_LIBS})
target_link_libraries(chanfast PRIVATE ${FFTW_LIBRARIES})
target_link_libraries(chanfast PRIVATE chanfast_lib)
target_link_libraries(chanfast PRIVATE m)
target_compile_features(chanfast PRIVATE cxx_std_17)
