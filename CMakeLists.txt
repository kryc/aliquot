cmake_minimum_required(VERSION 3.28)

# specify clang
SET(CMAKE_C_COMPILER clang)
SET(CMAKE_CXX_COMPILER clang++)

project(aliquot VERSION 0.1.0 LANGUAGES CXX)

enable_testing()

# Set build type default Release
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

add_subdirectory(./test/)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_COMPILER clang++ CACHE STRING "C++ compiler" FORCE)
# set(CMAKE_CXX_COMPILER /usr/bin/clang++)
set(CMAKE_C_COMPILER clang CACHE STRING "C compiler" FORCE)
# set(CMAKE_C_COMPILER /usr/bin/clang)

set(CMAKE_CXX_FLAGS "-Wall -Werror -Wno-unused-function -Wno-deprecated-declarations -Wno-c23-extensions")
if(!APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunsafe-buffer-usage -Wunsafe-buffer-usage-in-container")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-exceptions")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops -ffast-math -fomit-frame-pointer -fno-rtti")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -O0")
set(CMAKE_EXE_LINKER_FLAGS "-Wc23-extensions")
if(APPLE)
    execute_process(COMMAND brew --prefix OUTPUT_VARIABLE HOMEBREW_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
    message("Homebrew prefix: ${HOMEBREW_PREFIX}")
    include_directories(${HOMEBREW_PREFIX}/include)
    link_directories(${HOMEBREW_PREFIX}/lib)
endif()

set(ALIQUOT_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/aliquot.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/primefactors.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/primes.cpp
)
add_executable(aliquot ${ALIQUOT_SOURCES})
target_include_directories(aliquot
                            PRIVATE
                                ${CMAKE_CURRENT_SOURCE_DIR}/src
                        )
target_link_libraries(aliquot PRIVATE gmp gmpxx)

set(PRIMEGEN_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/primegen.cpp
)
add_executable(primegen ${PRIMEGEN_SOURCES})
target_include_directories(primegen
                            PRIVATE
                                ${CMAKE_CURRENT_SOURCE_DIR}/src
                        )
target_link_libraries(primegen PRIVATE gmp gmpxx)

set(FACTORGEN_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/factorgen.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/primefactors.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/primes.cpp
)
add_executable(factorgen ${FACTORGEN_SOURCES})
target_include_directories(factorgen
                            PRIVATE
                                ${CMAKE_CURRENT_SOURCE_DIR}/src
                        )
target_link_libraries(factorgen PRIVATE gmp gmpxx)

set(CACHECHECK_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cachecheck.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/primefactors.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/primes.cpp
)
add_executable(cachecheck ${CACHECHECK_SOURCES})
target_include_directories(cachecheck
                            PRIVATE
                                ${CMAKE_CURRENT_SOURCE_DIR}/src
                        )
target_link_libraries(cachecheck PRIVATE gmp gmpxx)

set(CACHESORT_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cachesort.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/primefactors.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/primes.cpp
)
add_executable(cachesort ${CACHESORT_SOURCES})
target_include_directories(cachesort
                            PRIVATE
                                ${CMAKE_CURRENT_SOURCE_DIR}/src
                        )
target_link_libraries(cachesort PRIVATE gmp gmpxx)