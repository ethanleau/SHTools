# CMakeList.txt : CMake project for SHTools, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
file(GLOB SOURCES
	 SHProjector.cpp
	 stb_image.h
	 stb_image_write.h
	 SH.h
	 SHEval.h
	 )

add_executable (SHTools ${SOURCES})

# TODO: Add tests and install targets if needed.
