set(CPLEX_ROOT_DIR "/opt/ibm/ILOG/CPLEX_Studio129/cplex/")
find_path(CPLEX_INCLUDE_DIR cplex.h HINTS "${CPLEX_ROOT_DIR}/include/ilcplex")
find_library(CPLEX_LIBRARY libcplex.a HINTS "${CPLEX_ROOT_DIR}/lib/x86-64_linux/static_pic")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CPLEX DEFAULT_MSG CPLEX_LIBRARY CPLEX_INCLUDE_DIR)

if(CPLEX_FOUND)
    set(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR})
    set(CPLEX_LIBRARIES ${CPLEX_LIBRARY})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(CPLEX_LIBRARIES "${CPLEX_LIBRARIES};m;pthread;dl")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(CPLEX_FOUND)

mark_as_advanced(CPLEX_LIBRARY CPLEX_INCLUDE_DIR)
