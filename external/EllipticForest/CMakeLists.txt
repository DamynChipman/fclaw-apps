# include(FetchContent)

option(ELLIPTICFOREST_PATH "Path to pre-installed EllipticForest" "")
if (ELLIPTICFOREST_PATH)

    if(EXISTS ${ELLIPTICFOREST_PATH}/include/EllipticForestApp.hpp)
        set(ELLIPTICFOREST_INCLUDE ${ELLIPTICFOREST_PATH}/include CACHE INTERNAL "EllipticForest include directory")
    else()
        message(FATAL_ERROR "Could not find EllipticForestApp.hpp")
    endif()

    if (EXISTS ${ELLIPTICFOREST_PATH}/lib/libEllipticForest.a)
        set(ELLIPTICFOREST_LIB ${ELLIPTICFOREST_PATH}/lib CACHE INTERNAL "EllipticForest library directory")
    else()
        message(FATAL_ERROR "Could not find libEllipticForest.a")
    endif()

    set(ELLIPTICFOREST_FOUND TRUE CACHE INTERNAL "Flag for if EllipticForest was found")

    message(STATUS "Found external EllipticForest:")
    message(STATUS "  ELLIPTICFOREST_INCLUDE = " ${ELLIPTICFOREST_INCLUDE})
    message(STATUS "  ELLIPTICFOREST_LIB = " ${ELLIPTICFOREST_LIB})

else()

    message(FATAL_ERROR "External EllipticForest not found or not set...")

endif()