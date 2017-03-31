set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFiniteCellApplication ****")

######### add include directory here
include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications )
include_directories( ${CMAKE_SOURCE_DIR}/applications/mkl_solvers_application/external_includes )
#######################################################################

######### add compilation units here
set( KRATOS_FINITE_CELL_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/dummy_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/dummy_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/element_wrapper_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/finite_cell_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_algebra_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/finite_cell_python_application.cpp
)
#######################################################################

## define library Kratos which defines the basic python interface
add_library(KratosFiniteCellApplication SHARED ${KRATOS_FINITE_CELL_APPLICATION_SOURCES})
target_link_libraries(KratosFiniteCellApplication KratosCore)
target_link_libraries(KratosFiniteCellApplication ${LAPACK_LIBRARIES})
set_target_properties(KratosFiniteCellApplication PROPERTIES PREFIX "")
install(TARGETS KratosFiniteCellApplication DESTINATION libs )
#######################################################################

## changing the .dll suffix to .pyd for Windows compatibility
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosFiniteCellApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#######################################################################

## install application's python script if required
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)
#######################################################################

## Add to the KratosMultiphysics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/FiniteCellApplication.py" DESTINATION KratosMultiphysics )
#######################################################################
