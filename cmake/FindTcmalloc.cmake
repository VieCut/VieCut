# - Find Tcmalloc
# Find the native Tcmalloc includes and library
#
#  Tcmalloc_INCLUDE_DIR - where to find Tcmalloc.h, etc.
#  Tcmalloc_LIBRARIES   - List of libraries when using Tcmalloc.
#  Tcmalloc_FOUND       - True if Tcmalloc found.

find_path(Tcmalloc_INCLUDE_DIR google/tcmalloc.h NO_DEFAULT_PATH PATHS
  ${HT_DEPENDENCY_INCLUDE_DIR}
  /usr/include
  /opt/local/include
  /usr/local/include
)

message(STATUS ${Tcmalloc_INCLUDE_DIR})

set(Tcmalloc_NAMES "")

message(STATUS ${Tcmalloc_NAMES})

find_library(Tcmalloc_LIBRARY
  NAMES tcmalloc libtcmalloc libtcmalloc_and_profiler libtcmalloc_minimal
  PATHS ${HT_DEPENDENCY_LIB_DIR} /lib /usr/lib /usr/local/lib /opt/local/lib /usr/lib/x86_84-linux-gnu
)

message(STATUS ${Tcmalloc_LIBRARY})

if (Tcmalloc_INCLUDE_DIR AND Tcmalloc_LIBRARY)
  set(Tcmalloc_FOUND TRUE)
  set( Tcmalloc_LIBRARIES ${Tcmalloc_LIBRARY} )
else ()
  set(Tcmalloc_FOUND FALSE)
  set( Tcmalloc_LIBRARIES )
endif ()

if (Tcmalloc_FOUND)
  message(STATUS "Found Tcmalloc: ${Tcmalloc_LIBRARY}")
else ()
  message(STATUS "Not Found Tcmalloc: ${Tcmalloc_LIBRARY}")
  if (Tcmalloc_FIND_REQUIRED)
    message(STATUS "Looked for Tcmalloc libraries named ${Tcmalloc_NAMES}.")
    message(FATAL_ERROR "Could NOT find Tcmalloc library")
  endif ()
endif ()

mark_as_advanced(
  Tcmalloc_LIBRARY
  Tcmalloc_INCLUDE_DIR
  )

