cmake_minimum_required(VERSION 3.06)

project(CG_chain)

option(G "Set to ON to compile with optimisations and debug symbols" OFF)

# The default behaviour is to compile with optimisations
set(CMAKE_BUILD_TYPE Release)

if(G)
	set(CMAKE_BUILD_TYPE RelWithDebInfo)
	message(STATUS "Compiling in Release+Debug symbols mode")
endif(G)

add_executable(Veff_chain Veff_chain.cpp utils/strings.cpp utils/Input.cpp)

add_executable(bead_chain bead_chain.cpp CG/LookupTable.cpp utils/strings.cpp utils/Input.cpp)
