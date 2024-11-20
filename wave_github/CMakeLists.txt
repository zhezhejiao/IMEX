
CMAKE_MINIMUM_REQUIRED(VERSION 3.10)


SET(TARGET "wave")
PROJECT(${TARGET})

SET(TARGET_SOURCES main/Main.cc)


FILE(GLOB TARGET_SRC  
    "local_integration/*.cpp"
    "exact_solution/*.cpp"
    "intial_value/*.cpp"
    "time_discretization/*.cpp"
    "InverseMatrix/*.cpp"
    "main/*.cpp"
)

FILE(GLOB TARGET_INC  
    "local_integration/*.h"
    "exact_solution/*.h"
    "intial_value/*.h"
    "time_discretization/*.h"
    "InverseMatrix/*.h"
    "main/*.h"
)


SET(TARGET_SRC 
    ${TARGET_SOURCES} 
    ${TARGET_SRC} 
    ${TARGET_INC}
)


FIND_PACKAGE(deal.II 8.3 QUIET
    HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
)

IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
  )
ENDIF()


DEAL_II_INITIALIZE_CACHED_VARIABLES()


DEAL_II_INVOKE_AUTOPILOT()
