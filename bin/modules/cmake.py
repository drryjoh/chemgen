import yaml
import os

def generate_cmake_file(config, output_path):
    """Generate a CMakeLists.txt file based on the configuration."""
    cmake_content = f"""\
cmake_minimum_required(VERSION 3.15)
project(ChemGen)  # Replace with your project name

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Compiler options
set(CMAKE_CXX_FLAGS "-{config['build'].get('chemgen_optimized', 'O2')}")

# Source files
set(SOURCE_FILES
    ${{CMAKE_SOURCE_DIR}}/src/main.cpp   # Add your source files here
)

# Set the output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${{CMAKE_BINARY_DIR}}/bin)

# Add the executable
add_executable(output_program ${{SOURCE_FILES}})

# Link libraries based on configuration
"""

    # Add TBB linking if specified
    if config['build'].get('chemgen_smp') == 'TBB':
        cmake_content += """\
find_package(TBB REQUIRED)
target_link_libraries(output_program PRIVATE TBB::tbb)
"""

    # Enable debugging flags if needed
    cmake_content += "\n# Enable debugging flags\nset(CMAKE_BUILD_TYPE Debug)"

    # Write content to CMakeLists.txt
    with open(output_path/"CMakeLists.txt", 'w') as file:
        file.write(cmake_content)

# Load configuration and generate CMakeLists.txt
