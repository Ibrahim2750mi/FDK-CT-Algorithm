cmake_minimum_required(VERSION 3.10)
project(FDK_Reconstruction)

set(CMAKE_CXX_STANDARD 11)

find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED fftw3)

find_package(OpenCV REQUIRED)

add_executable(fdk_recon
        src/main.cpp
        src/fdk.cpp
        # src/phantom.cpp
)

# Include directories
target_include_directories(fdk_recon PRIVATE
        ${FFTW_INCLUDE_DIRS}
        ${OpenCV_INCLUDE_DIRS}
)

# Link libraries
target_link_libraries(fdk_recon
        ${FFTW_LIBRARIES}
        ${OpenCV_LIBS}
)
