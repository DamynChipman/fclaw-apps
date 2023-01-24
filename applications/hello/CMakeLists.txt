add_executable(hello main.cpp)
target_include_directories(hello PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${EXTERNAL_PKGS_INCLUDE})
target_link_directories(hello PUBLIC ${EXTERNAL_PKGS_LIB})
target_link_libraries(hello PUBLIC z sc p4est petsc forestclaw EllipticForest)