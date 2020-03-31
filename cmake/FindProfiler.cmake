# - Find Tcmalloc Profiler
# Find the native Tcmalloc Profiler includes and library
#
#  Profiler_INCLUDE_DIR - where to find Profiler.h, etc.
#  Profiler_LIBRARIES   - List of libraries when using Profiler.
#  Profiler_FOUND       - True if Profiler found.

# Adapted from https://github.com/COMBINE-lab/quark/blob/master/cmake/Modules/FindTcmalloc.cmake

find_path(Profiler_INCLUDE_DIR google/tcmalloc.h NO_DEFAULT_PATH PATHS
  $ENV{HOME}/.local/include
  /usr/include
  /opt/local/include
  /usr/local/include
)

set(Profiler_NAMES "")

find_library(Profiler_LIBRARY
  NAMES libtcmalloc_and_profiler tcmalloc_and_profiler NO_DEFAULT_PATH
  PATHS $ENV{HOME}/.local/lib
        /lib
        /usr/lib
        /usr/local/lib
        /opt/local/lib
        /usr/lib/x86_84-linux-gnu
)

if (Profiler_INCLUDE_DIR AND Profiler_LIBRARY)
  set(Profiler_FOUND TRUE)
  set( Profiler_LIBRARIES ${Profiler_LIBRARY} )
else ()
  set(Profiler_FOUND FALSE)
  set( Profiler_LIBRARIES )
endif ()

if (Profiler_FOUND)
  message(STATUS "Found Profiler: ${Profiler_LIBRARY} ${Profiler_INCLUDE_DIR}")
else ()
  message(STATUS "Not Found Profiler: ${Profiler_LIBRARY}")
  if (Profiler_FIND_REQUIRED)
    message(STATUS "Looked for Profiler libraries named ${Profiler_NAMES}.")
    message(FATAL_ERROR "Could NOT find Profiler library")
  endif ()
endif ()

mark_as_advanced(
  Profiler_LIBRARY
  Profiler_INCLUDE_DIR
  )
