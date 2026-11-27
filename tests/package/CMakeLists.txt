# Projet consommateur, lancé par la CI : vérifie que wsl_libraries s'utilise
# comme le README l'annonce. MODE = package (find_package après installation)
# ou subdirectory (FetchContent / add_subdirectory sur le dossier source).
cmake_minimum_required(VERSION 3.28)
project(wsl_libraries_consumer LANGUAGES CXX)

# Volontairement C++20 : lier wsl::mdarray doit faire passer en C++23.
set(CMAKE_CXX_STANDARD 20)

set(MODE "package" CACHE STRING "package ou subdirectory")
if(MODE STREQUAL "package")
    find_package(wsl_libraries 0.1 REQUIRED)
else()
    include(FetchContent)
    FetchContent_Declare(wsl_libraries SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../..)
    FetchContent_MakeAvailable(wsl_libraries)
    if(TARGET test_permutation)
        message(FATAL_ERROR "Les tests de wsl_libraries ne doivent pas être construits chez un consommateur")
    endif()
endif()

add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE wsl::mdarray wsl::permutation wsl::random wsl::special)

enable_testing()
add_test(NAME consumer COMMAND consumer)
