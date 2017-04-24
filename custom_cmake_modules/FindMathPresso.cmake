#
# Looks for MathPresso packages
#

IF(MATHPRESSO_DIR)
    SET(MATHPRESSO_INC_DIRS ${MATHPRESSO_DIR}/include)
    SET(MATHPRESSO_LIB_DIRS ${MATHPRESSO_DIR}/lib)

    INCLUDE_DIRECTORIES(${MATHPRESSO_INC_DIRS})

    FIND_LIBRARY(AUX mathpresso ${MATHPRESSO_LIB_DIRS} NO_DEFAULT_PATH)
    
    IF(AUX)
        SET(MATHPRESSO_LIBRARIES ${AUX})
        SET(MATHPRESSO_FOUND TRUE)
    ELSE(AUX)
        SET(MATHPRESSO_FOUND FALSE)
    ENDIF(AUX)
ELSE(MATHPRESSO_DIR)
    SET(MATHPRESSO_FOUND FALSE)
ENDIF(MATHPRESSO_DIR)

