cmake_minimum_required(VERSION 3.20)
project(intersection_gen LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_compile_definitions(PROJECT_ROOT_DIR=\"${PROJECT_SOURCE_DIR}\")

set(DEPENDS_EXISTS "FALSE")
set(DEPENDS_DIR "${CMAKE_SOURCE_DIR}/depends")
if(EXISTS "${DEPENDS_DIR}")
    set(DEPENDS_EXISTS "TRUE")
    message(STATUS "Found [DEPENDS_DIR=${DEPENDS_DIR}], using it directly.")
    ## 添加第三方依赖库目录
    include_directories(${DEPENDS_DIR}/nlohmann/)
    include_directories(${DEPENDS_DIR}/eigen/include/eigen3)
#    include_directories(${DEPENDS_DIR}/clipper2/include)
    #include_directories(${DEPENDS_DIR}/geos/include)

    link_directories(${DEPENDS_DIR}/eigen/lib)
#    link_directories(${DEPENDS_DIR}/clipper2/lib)
    #link_directories(${DEPENDS_DIR}/geos/lib)
else()
    set(DEPENDS_EXISTS "FALSE")
    message(STATUS "Not found [DEPENDS_DIR=${DEPENDS_DIR}], using <FetchContent> fetch libriaries.")
    include(FetchContent)
#    ## nlohmann_json
#    message(STATUS "Start fetch libriary: nlohmann_json ...")
#    FetchContent_Declare(nlohmann_json
#            GIT_REPOSITORY https://github.com/nlohmann/json.git
#            GIT_TAG        v3.11.3)
#    FetchContent_MakeAvailable(nlohmann_json)
    ## eigen
    message(STATUS "Start fetch libriary: eigen ...")
    FetchContent_Declare(eigen
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
            GIT_TAG        3.4.0
            GIT_SHALLOW    TRUE)
    set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING   OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(eigen)
#    ## Clipper2
#    message(STATUS "Start fetch libriary: Clipper2 ...")
#    option(CLIPPER2_EXAMPLES "Build examples" OFF)
#    option(CLIPPER2_TESTS "Build tests" OFF)
#    option(USE_EXTERNAL_GTEST "Use system-wide installed GoogleTest" OFF)
#    option(USE_EXTERNAL_GBENCHMARK "Use the googlebenchmark" OFF)
#    FetchContent_Declare(
#            Clipper2
#            GIT_REPOSITORY https://github.com/AngusJohnson/Clipper2.git
#            GIT_TAG        main
#            SOURCE_SUBDIR  CPP  # change directory to CMakeLists.txt
#            GIT_SHALLOW    TRUE)  # ver: main | Clipper2_1.4.0
#    FetchContent_MakeAvailable(Clipper2)
endif()

add_subdirectory(src)
add_subdirectory(tests)
