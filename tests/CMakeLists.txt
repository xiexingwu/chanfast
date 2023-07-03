# Google Test
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.13.0
)
FetchContent_MakeAvailable(googletest)

add_executable(tests
  test_decomp.cpp
)
target_include_directories(tests PRIVATE ${chanfast_INCLUDES})
target_include_directories(tests PRIVATE ${FFTW_INCLUDE_DIRS})
target_link_libraries(tests PRIVATE ${chanfast_LIBS})
target_link_libraries(tests PRIVATE ${FFTW_LIBRARIES})
target_link_libraries(tests PRIVATE chanfast_lib)
target_link_libraries(tests PRIVATE gtest_main)
target_compile_features(tests PRIVATE cxx_std_17)
