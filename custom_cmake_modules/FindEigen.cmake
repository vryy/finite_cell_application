#
# Looks for EIGEN package
#

message("EIGEN_DIR:" ${EIGEN_DIR})
IF(EIGEN_DIR)
    SET(EIGEN_INC_DIRS ${EIGEN_DIR})

    INCLUDE_DIRECTORIES(${EIGEN_INC_DIRS})

    SET(EIGEN_FOUND TRUE)
ELSE()
    SET(EIGEN_FOUND FALSE)
ENDIF()
