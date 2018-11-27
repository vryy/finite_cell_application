#
# Looks for SPECTRA package
#

IF(${EIGEN_FOUND} MATCHES TRUE)
    message("SPECTRA_DIR:" ${SPECTRA_DIR})
    IF(SPECTRA_DIR)
        SET(SPECTRA_INC_DIRS ${SPECTRA_DIR}/include)
        INCLUDE_DIRECTORIES(${SPECTRA_INC_DIRS})
        SET(SPECTRA_FOUND TRUE)
    ELSE()
        SET(SPECTRA_FOUND FALSE)
    ENDIF()
ELSE()
    message("Eigen is not found. Spectra dependencies are not met")
    SET(SPECTRA_FOUND FALSE)
ENDIF()

