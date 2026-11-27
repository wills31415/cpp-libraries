# Projet consommateur, lancé par la CI : vérifie que le paquet wsl s'utilise
# comme le README l'annonce. MODE = package (find_package après installation)
# ou subdirectory (FetchContent / add_subdirectory sur le dossier source).
cmake_minimum_required(VERSION 3.28)
project(wsl_consumer LANGUAGES CXX)

# Volontairement C++20 : lier wsl::mdarray doit faire passer en C++23.
set(CMAKE_CXX_STANDARD 20)

set(MODE "package" CACHE STRING "package ou subdirectory")
if(MODE STREQUAL "package")
    find_package(wsl 0.1 REQUIRED)
else()
    include(FetchContent)
    FetchContent_Declare(wsl SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../..)
    FetchContent_MakeAvailable(wsl)
    if(TARGET permutation)
        message(FATAL_ERROR "Les tests de wsl ne doivent pas être construits chez un consommateur")
    endif()
endif()

add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE wsl::mdarray wsl::permutation wsl::random wsl::special)
# La cible de tout le dépôt.
add_executable(consumer_all main.cpp)
target_link_libraries(consumer_all PRIVATE wsl::wsl)

enable_testing()
add_test(NAME consumer COMMAND consumer)
add_test(NAME consumer_all COMMAND consumer_all)
