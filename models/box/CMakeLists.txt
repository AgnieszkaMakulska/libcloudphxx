cmake_minimum_required(VERSION 3.26)
project(box)
enable_testing()
set(CMAKE_CXX_STANDARD 17)
add_executable(box box.cpp)


find_package(libcloudph++)
if (NOT libcloudph++_FOUND)
    message(FATAL_ERROR "libclouphxx++ not found")
endif()

target_compile_features(box PUBLIC cxx_std_20)

target_link_libraries(box PRIVATE clphxx::cloudphxx_lgrngn)

find_package(Boost COMPONENTS program_options REQUIRED)
target_link_libraries(box PRIVATE ${Boost_LIBRARIES})
target_include_directories(box PUBLIC ${Boost_INCLUDE_DIRS})


find_library(netcdfcxx_LIBRARIES netcdf_c++4 REQUIRED)
target_link_libraries(box PRIVATE ${netcdfcxx_LIBRARIES})
find_path(netcdfxx_INCLUDE_DIRS netcdf)
target_include_directories(box PUBLIC ${netcdfxx_INCLUDE_DIRS})

find_package(Thrust REQUIRED)
thrust_create_target(Thrust)
target_link_libraries(box PRIVATE Thrust)