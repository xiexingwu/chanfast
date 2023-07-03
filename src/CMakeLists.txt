set(chanfastlib_SOURCES
  mpi_plan.cpp
  decomp.cpp
  tensor.cpp
)

add_library(chanfast_lib STATIC ${chanfastlib_SOURCES})
target_include_directories(chanfast_lib PRIVATE ${chanfast_INCLUDES})
target_include_directories(chanfast_lib PRIVATE ${FFTW_INCLUDE_DIRS})
target_compile_features(chanfast_lib PRIVATE cxx_std_17)
