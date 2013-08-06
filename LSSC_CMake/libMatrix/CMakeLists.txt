CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(libMatrix)

SET(SRC LibMatrix.h LibMatrix.cpp ClassMatrix.h ClassMatrix.cpp)

ADD_LIBRARY(MATRIX ${SRC})

INCLUDE_DIRECTORIES(..)		
	
INSTALL(FILES CMakeLists.txt ${SRC} DESTINATION src/libMatrix)
