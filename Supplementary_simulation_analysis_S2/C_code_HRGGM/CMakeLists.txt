cmake_minimum_required(VERSION 3.16)

project(C_GRN LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(C_GRN
    Main.cpp
    Matrix_io.cpp
    Matrix_compute.cpp
    LHG_update.cpp
    Class_define.cpp
    Parameter_select_parallel.cpp
    Newton_method_offdiag_l.cpp
    Newton_method_offdiag_lhg.cpp
)

target_include_directories(C_GRN PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Eigen
)

find_package(OpenMP REQUIRED)
target_link_libraries(C_GRN PRIVATE OpenMP::OpenMP_CXX)

if (MSVC)
    target_compile_options(C_GRN PRIVATE /O2)
else()
    target_compile_options(C_GRN PRIVATE -O2 -Wall -Wextra)
endif()