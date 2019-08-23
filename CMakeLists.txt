set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFiniteCellApplication ****")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/custom_cmake_modules")

######### add include directory here
include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications )
include_directories( ${CMAKE_SOURCE_DIR}/applications/mkl_solvers_application/external_includes )
include_directories( ${CMAKE_SOURCE_DIR}/applications/brep_application )

if(${ISOGEOMETRIC_APPLICATION} MATCHES ON)
    include_directories( ${CMAKE_SOURCE_DIR}/applications/isogeometric_application )
    ADD_DEFINITIONS(-DENABLE_FINITE_CELL_ISOGEOMETRIC)
endif()

if(EXISTS ${CMAKE_SOURCE_DIR}/applications/mkl_solvers_application/external_includes)
    ADD_DEFINITIONS(-DENABLE_FINITE_CELL_BOOST_BINDINGS)
endif()

#######################################################################

FIND_PACKAGE(NNLS REQUIRED)
IF(${NNLS_FOUND} MATCHES TRUE)
    MESSAGE("nnls is found at " ${NNLS_DIR})
    ADD_DEFINITIONS(-DFINITE_CELL_APPLICATION_USE_NNLS)
ELSE()
    MESSAGE("nnls is not found")
ENDIF()

######### add compilation units here
set( KRATOS_FINITE_CELL_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/dummy_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/dummy_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/element_wrapper_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/quad_tree_node.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/ghost_penalty_utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/finite_cell_mesh_utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/finite_cell_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_conditions_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_quadrature_utility_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_div_free_basis_utility_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_finite_cell_mesh_utility_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_finite_cell_auxiliary_utility_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_moment_fitting_utility_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_ghost_penalty_utility_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_utility_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/finite_cell_python_application.cpp
)
#######################################################################

## define library Kratos which defines the basic python interface
add_library(KratosFiniteCellApplication SHARED ${KRATOS_FINITE_CELL_APPLICATION_SOURCES})
target_link_libraries(KratosFiniteCellApplication KratosCore)
target_link_libraries(KratosFiniteCellApplication KratosBRepApplication)
target_link_libraries(KratosFiniteCellApplication ${LAPACK_LIBRARIES})
IF(${NNLS_FOUND} MATCHES TRUE)
    target_link_libraries(KratosFiniteCellApplication ${NNLS_LIBRARIES})
ENDIF()
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
