cmake_minimum_required(VERSION 3.5)

# name of the project
project(ttk_cpp)
set(Boost_INCLUDE_DIR Path_to_Boost_1_73)

if (APPLE)
  execute_process(COMMAND brew --prefix libomp
                OUTPUT_VARIABLE OpenMP_HOME
                OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "OpenMP Root : ${OpenMP_HOME}")

  set(OpenMP_C_LIB_NAMES "omp")
  set(OpenMP_CXX_LIB_NAMES "omp")
  set(OpenMP_omp_LIBRARY "${OpenMP_HOME}/lib/")
  set(OpenMP_CXX_FLAGS "-Xpreprocessor -I${OpenMP_HOME}/include" CACHE STRING "" FORCE)
  set(OpenMP_C_FLAGS "-I${OpenMP_HOME}/include" CACHE STRING "" FORCE)
endif (APPLE)

set(TTKBase_DIR Path_to_ttkBase)
set(CMAKE_CXX_STANDARD 11)
unset(TTK_ENABLE_OPENMP CACHE)
find_package(TTKBase REQUIRED)
find_package(OpenMP REQUIRED)


set(TARGET ttk_cpp)
add_library(${TARGET}
  STATIC
    ttk_cpp.cpp
)

target_link_libraries(${TARGET}
  PRIVATE
    ttk::base::baseAll
    )


install (TARGETS ${TARGET} DESTINATION .)
