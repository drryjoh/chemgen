import yaml
import os

def generate_cmake_file(config, output_path, third_parties):
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
    ${{CMAKE_SOURCE_DIR}}/src/chemgen.cpp   # Add your source files here
)

# Set the output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${{CMAKE_BINARY_DIR}}/bin)

# Add the executable
add_executable(chemgen ${{SOURCE_FILES}})

# Link libraries based on configuration
"""

    # Add TBB linking if specified
    [use_third_parties, third_party_path, libraries] = third_parties


    if use_third_parties:
        cmake_content+="""\
# Include TBB as an external project
include(ExternalProject)
        """
        for library in libraries:
            cmake_content += f"""\

# Configure the {library.upper()} build
ExternalProject_Add({library}
    SOURCE_DIR {third_party_path}/{library}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${{CMAKE_BINARY_DIR}}/{library}_install
        -D{library.upper()}_TEST=OFF
        -D{library.upper()}_EXAMPLES=OFF
        -D{library.upper()}_ENABLE_IPO=OFF
    BUILD_IN_SOURCE 1
    UPDATE_DISCONNECTED 1
    BUILD_COMMAND ${{CMAKE_COMMAND}} --build . --config Release
    INSTALL_COMMAND ${{CMAKE_COMMAND}} --build . --target install
)

# Set the TBB include and library paths
set({library.upper()}_INCLUDE_DIR ${{CMAKE_BINARY_DIR}}/{library}_install/include)
set({library.upper()}_LIBRARY_DIR ${{CMAKE_BINARY_DIR}}/{library}_install/lib)

# Include directories for the project
include_directories(${{{library.upper()}_INCLUDE_DIR}})

# Link TBB libraries
target_link_libraries(chemgen PRIVATE ${{{library.upper()}_LIBRARY_DIR}}/lib{library}.so)

# Ensure TBB is built before chemgen
add_dependencies(chemgen {library})
"""

    # Enable debugging flags if needed
    cmake_content += "\n# Enable debugging flags\nset(CMAKE_BUILD_TYPE Release)"

    # Write content to CMakeLists.txt
    with open(output_path/"CMakeLists.txt", 'w') as file:
        file.write(cmake_content)

# Load configuration and generate CMakeLists.txt
