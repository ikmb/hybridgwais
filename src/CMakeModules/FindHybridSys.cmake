
set(HybridSys_INC_PATHS /usr/include /opt/biowulf/include hybridsys-api ../hybridsys-api)
set(HybridSys_LIB_PATHS /usr/lib)

find_library(HybridSys_LIBRARIES NAMES biowulf)
find_path(HybridSys_INCLUDE_DIRS biowulf.h PATHS ${ADMXRC3_INC_PATHS})

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(HybridSys DEFAULT_MSG HybridSys_LIBRARIES HybridSys_INCLUDE_DIRS)

