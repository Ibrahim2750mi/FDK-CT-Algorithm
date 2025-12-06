cmake_minimum_required(VERSION 3.10)
project(FDK_Reconstruction)

set(CMAKE_CXX_STANDARD 11)

find_package(PkgConfig REQUIRED)

find_package(OpenCV REQUIRED)

add_executable(fdk_recon
        src/main.cpp
        src/fdk.cpp
        # src/phantom.cpp
)

# Include directories
target_include_directories(fdk_recon PRIVATE
        ${OpenCV_INCLUDE_DIRS}
)

# Link libraries
target_link_libraries(fdk_recon
        ${OpenCV_LIBS}
)
