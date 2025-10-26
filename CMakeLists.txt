cmake_minimum_required(VERSION 3.30)

project(usearch LANGUAGES C CXX VERSION 12.0)

set(BUILD_SHARED_LIBS OFF)

include(FetchContent)

# configure zlib
set(ZLIB_BUILD_EXAMPLES OFF)

FetchContent_Declare(
    zlib
    GIT_REPOSITORY https://github.com/madler/zlib.git
    GIT_TAG 51b7f2abdade71cd9bb0e7a373ef2610ec6f9daf # version 1.3.1
)
FetchContent_MakeAvailable(zlib)

file(GLOB SRC_FILES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
)

# link with pthreads
find_package(Threads REQUIRED)

add_executable(usearch12)
target_sources(usearch12 PRIVATE ${SRC_FILES})
target_link_libraries(usearch12 PRIVATE zlibstatic Threads::Threads)

target_compile_features(usearch12 PUBLIC cxx_std_11)

if(MSVC)
    target_compile_definitions(usearch12 PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()
if(UNIX)
    target_compile_options(usearch12 PRIVATE -ffast-math)
endif()
