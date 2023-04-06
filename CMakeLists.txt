project(GI)
CMAKE_MINIMUM_REQUIRED(VERSION 3.16)


set(GMP_INCLUDE_DIR 3rdparty/gmp/include)
set(GMP_LIBRARIES 3rdparty/gmp/lib)

set(MPFR_INCLUDE_DIR 3rdparty/gmp/include)
set(MPFR_LIBRARIES 3rdparty/gmp/lib)

set(CGAL_DIR 3rdparty/cgal)
find_package(CGAL REQUIRED)

set(OpenCV_DIR 3rdparty/opencv/build)
find_package(OpenCV REQUIRED)


file(GLOB SRCLIST ./src/*.cpp ./src*.hpp)
if (MSVC)
  add_compile_options(/bigobj)
else ()
  add_compile_options(-Wa,-mbig-obj)
endif ()
add_definitions(-D_USE_MATH_DEFINES)

include_directories(include)
include_directories(3rdparty/eigen3)
include_directories(3rdparty/libigl/include)
include_directories(3rdparty/boost)

include_directories(${GMP_INCLUDE_DIR})
include_directories(${MPFR_INCLUDE_DIR})
include_directories(${CGAL_INCLUDE_DIRS})
include_directories(${OpenCV_INCLUDE_DIRS})

link_directories(C:/Boost/lib)

add_executable(GI ${SRCLIST})

target_link_libraries(GI libboost_filesystem-vc142-mt-gd-x64-1_81.lib 
                         ${OpenCV_LIBRARIES})



