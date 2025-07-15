cmake_minimum_required(VERSION 3.30.0)

project(eigenfem1 VERSION 0.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(SOURCE_FILES
    ${PROJECT_SOURCE_DIR}/src/eigenfem/main.cpp
    ${PROJECT_SOURCE_DIR}/src/eigenfem/element.h
    ${PROJECT_SOURCE_DIR}/src/eigenfem/element.cpp
    ${PROJECT_SOURCE_DIR}/src/eigenfem/materials.h
    ${PROJECT_SOURCE_DIR}/src/eigenfem/materials.cpp
    ${PROJECT_SOURCE_DIR}/src/eigenfem/mesh.h
    ${PROJECT_SOURCE_DIR}/src/eigenfem/mesh.cpp
    ${PROJECT_SOURCE_DIR}/src/eigenfem/model.h
    ${PROJECT_SOURCE_DIR}/src/eigenfem/model.cpp
    ${PROJECT_SOURCE_DIR}/src/eigenfem/solvers.h
    ${PROJECT_SOURCE_DIR}/src/eigenfem/solvers.cpp
    ${PROJECT_SOURCE_DIR}/src/eigenfem/utils.cpp
    ${PROJECT_SOURCE_DIR}/src/eigenfem/utils.h
    ${PROJECT_SOURCE_DIR}/src/eigenfem/io.h
    ${PROJECT_SOURCE_DIR}/src/eigenfem/io.cpp
    )

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

add_library(eigen INTERFACE)
add_library(eigen2 INTERFACE)
add_library(spectra INTERFACE)

target_include_directories(eigen INTERFACE ${PROJECT_SOURCE_DIR}/third-party/eigen-3.4.0/Eigen)
target_include_directories(eigen2 INTERFACE ${PROJECT_SOURCE_DIR}/third-party/eigen-3.4.0/)
target_include_directories(spectra INTERFACE ${PROJECT_SOURCE_DIR}/third-party/spectra-1.1.0/include)
target_link_libraries(${PROJECT_NAME} PRIVATE eigen)
target_link_libraries(${PROJECT_NAME} PRIVATE eigen2)
target_link_libraries(${PROJECT_NAME} PRIVATE spectra)