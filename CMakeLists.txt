# @file HdivDarcy/CMakeList.txt  -- First file to read for CMake
project (HdivDarcy)

cmake_minimum_required(VERSION 3.11.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

IF (WIN32)
   find_package( PZ PATHS "C:/Arquivos de Programas/PZ" REQUIRED NO_DEFAULT_PATH)
ELSE()
   find_package( PZ REQUIRED )
ENDIF (WIN32)

include_directories(${PZ_INCLUDE_DIRS})

IF(APPLE)
   # Mac OS X specific code
   add_definitions(-DMACOSX)
   set(MACOSX MACOSX)
   FIND_LIBRARY(ACCELERATE_LIB Accelerate)
   link_libraries(${ACCELERATE_LIB})
ENDIF(APPLE)


# Location to search files or libraries
set(SEARCH_DIRS ENV PATH /usr/local/opt/ /usr/local/opt/qt5/lib/cmake/ /usr/include /usr/local/include /usr/local/include/boost /usr/local/include/log4cxx /usr/local/include/tbb /usr/lib /usr/lib64 /usr/local/lib ./. ./.. ./../libs ./../libs/fad ./../libs/metis ./../libs/metis/include ./../libs/metis/lib ./../externallibs ./../externallibs/fad ./../externallibs/pthread ./../externallibs/metis ./../externallibs/metis/include ./../externallibs/metis/lib ./..metis ./../metis/include ./../metis/lib ./externallibs/lib ./externallibs/include ./../externallibs/lib ./../externallibs/include)

if(WIN32)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")
   #add a define saying it is a VC compiler
   set(VC "#Will add a define VC on pz_config.h")
   add_definitions(-DVC)
   #define _USE_MATH_DEFINES for example: M_PI constant definitions.
   add_definitions(-D_USE_MATH_DEFINES)
   #disabling VC warnings
   add_definitions(-D_SCL_SECURE_NO_WARNINGS)
   add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	
   #define use of pthread static lib.
   add_definitions(-DPTW32_STATIC_LIB)
   #enabling /bigobj
   add_definitions("/bigobj")
 
   FOREACH(FLAG_TYPE EXE MODULE SHARED)
       STRING (REPLACE "INCREMENTAL:YES" "INCREMENTAL" FLAG_TMP "${CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG}")
       STRING (REPLACE "INCREMENTAL:NO" "INCREMENTAL" FLAG_TMP ${FLAG_TMP})
       STRING (REPLACE "INCREMENTAL" "INCREMENTAL:NO" FLAG_TMP ${FLAG_TMP})
       STRING (REPLACE "/EDITANDCONTINUE" "" FLAG_TMP ${FLAG_TMP})
       SET(CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG "${FLAG_TMP}" CACHE STRING "Overriding default debug ${FLAG_TYPE} linker flags." FORCE)
       MARK_AS_ADVANCED (CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG)
   ENDFOREACH ()
else()
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")#this flag is for preventing errors when casting from double to float
endif()

# Adding projec files
file(GLOB header *.h)
file(GLOB source *.cpp)

add_executable(HdivDarcy ${source} ${header})

target_link_libraries(HdivDarcy pz)

#copying input data folders
add_custom_command( TARGET HdivDarcy POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory  
	"${CMAKE_CURRENT_LIST_DIR}/input"
	$<TARGET_FILE_DIR:HdivDarcy>/input
)

add_custom_command( TARGET HdivDarcy POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory  
	"${CMAKE_CURRENT_LIST_DIR}/exp_data"
	$<TARGET_FILE_DIR:HdivDarcy>/exp_data
)