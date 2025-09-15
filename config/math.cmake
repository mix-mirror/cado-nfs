
set(CMAKE_REQUIRED_LIBRARIES m)

include(${PROJECT_SOURCE_DIR}/config/search_for_function.cmake)
search_for_function(log2 HAVE_LOG2)
search_for_function(exp2 HAVE_EXP2)

set(math_libs ${CMAKE_REQUIRED_LIBRARIES_EXTRA})
