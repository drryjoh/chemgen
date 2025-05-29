cmake_minimum_required(VERSION 3.15)
project(ChemGen)  # Replace with your project name

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Compiler options
set(CMAKE_CXX_FLAGS "-O2")

# Source files
set(SOURCE_FILES
    ${CMAKE_SOURCE_DIR}/src/chemgen.cpp   # Add your source files here
)

# Set the output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add the executable
add_executable(chemgen ${SOURCE_FILES})

# Link libraries based on configuration

# Enable debugging flags
set(CMAKE_BUILD_TYPE Release)