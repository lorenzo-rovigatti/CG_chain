cmake_minimum_required(VERSION 3.06)

# set the project name
project(CG_chain)

# The default behaviour is to compile with optimisations
set(CMAKE_BUILD_TYPE Release)

# add the executable
add_executable(Veff_chain Veff_chain.cpp utils/strings.cpp utils/Input.cpp)
