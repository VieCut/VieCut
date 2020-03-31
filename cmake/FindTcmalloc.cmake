# - Find Tcmalloc
# Find the native Tcmalloc includes and library
#
#  Tcmalloc_INCLUDE_DIR - where to find Tcmalloc.h, etc.
#  Tcmalloc_LIBRARIES   - List of libraries when using Tcmalloc.
#  Tcmalloc_FOUND       - True if Tcmalloc found.

# Adapted from https://github.com/COMBINE-lab/quark/blob/master/cmake/Modules/FindTcmalloc.cmake

find_path(Tcmalloc_INCLUDE_DIR google/tcmalloc.h NO_DEFAULT_PATH PATHS
  $ENV{HOME}/.local/include
  /usr/include
  /opt/local/include
  /usr/local/include
)

set(Tcmalloc_NAMES "")

find_library(Tcmalloc_LIBRARY
  NAMES tcmalloc libtcmalloc libtcmalloc_and_profiler libtcmalloc_minimal
  PATHS $ENV{HOME}/.local/lib
        /lib
        /usr/lib
        /usr/local/lib
        /opt/local/lib
        /usr/lib/x86_84-linux-gnu
)

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

