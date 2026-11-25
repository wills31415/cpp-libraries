cmake_minimum_required(VERSION 3.28)

project(wsl_libraries
        VERSION 0.1.0
        DESCRIPTION "Collection de bibliothèques C++ génériques, en en-têtes seuls"
        LANGUAGES CXX)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Tests et exemples : construits seulement quand ce dépôt est le projet principal.
# Un projet qui nous intègre par FetchContent ou add_subdirectory ne les voit pas.
option(WSL_LIBRARIES_BUILD_TESTS "Construire les tests" ${PROJECT_IS_TOP_LEVEL})
option(WSL_LIBRARIES_BUILD_EXAMPLES "Construire les exemples" ${PROJECT_IS_TOP_LEVEL})
option(WSL_LIBRARIES_SANITIZERS "Compiler tests et exemples avec ASan et UBSan" OFF)
option(WSL_LIBRARIES_WERROR "Traiter les avertissements des tests et exemples comme des erreurs (CI)" OFF)
# CI : fichier listant les chemins modifiés (un par ligne, « ALL » pour tout) ; seuls
# les tests et exemples des bibliothèques touchées, et de celles qui en dépendent, sont
# construits (cmake/affected.cmake). Vide : tout.
set(WSL_LIBRARIES_AFFECTED_BY "" CACHE FILEPATH "Fichier des chemins modifiés (CI)")

include(GNUInstallDirs)

# Standards sur lesquels les tests sont compilés. Une bibliothèque est testée
# sur chacun de ceux qu'elle prétend supporter (de son minimum à C++23). Le dépôt est
# en C++23 partout (docs/interfaces.md, décision 7) ; ajouter 20 ici si une
# bibliothèque redescend un jour à C++20, sur décision explicite.
set(WSL_LIBRARIES_STANDARDS 23)

# ---------------------------------------------------------------------------
# Les bibliothèques
# ---------------------------------------------------------------------------
#
#     wsl_add_library(<nom> STANDARD <20|23> [DEPENDS <nom>...])
#
# Déclare la bibliothèque en en-têtes seuls include/wsl/<nom>.hpp (plus, s'il
# existe, le dossier include/wsl/<nom>/ et ses sous-dossiers), sous la
# cible wsl::<nom>. STANDARD est le minimum exigé, transmis aux consommateurs
# par target_compile_features. DEPENDS nomme d'autres bibliothèques du dépôt.

set(WSL_LIBRARIES_ALL "")

function(wsl_add_library name)
    cmake_parse_arguments(PARSE_ARGV 1 ARG "" "STANDARD" "DEPENDS")
    if(NOT ARG_STANDARD)
        message(FATAL_ERROR "wsl_add_library(${name}) : STANDARD manquant")
    endif()

    set(target wsl_${name})
    add_library(${target} INTERFACE)
    add_library(wsl::${name} ALIAS ${target})
    set_target_properties(${target} PROPERTIES
        EXPORT_NAME ${name}
        WSL_STANDARD ${ARG_STANDARD})
    target_compile_features(${target} INTERFACE cxx_std_${ARG_STANDARD})
    file(GLOB_RECURSE headers CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/include/wsl/${name}/*.hpp)
    target_sources(${target} INTERFACE
        FILE_SET HEADERS
        BASE_DIRS ${PROJECT_SOURCE_DIR}/include
        FILES ${PROJECT_SOURCE_DIR}/include/wsl/${name}.hpp ${headers})
    foreach(dependency IN LISTS ARG_DEPENDS)
        target_link_libraries(${target} INTERFACE wsl_${dependency})
    endforeach()

    set_property(GLOBAL APPEND PROPERTY WSL_LIBRARY_NAMES ${name})
    set_property(GLOBAL PROPERTY WSL_DEPENDS_${name} ${ARG_DEPENDS})
    set(WSL_LIBRARIES_ALL ${WSL_LIBRARIES_ALL} ${target} PARENT_SCOPE)
endfunction()

define_property(TARGET PROPERTY WSL_STANDARD
    BRIEF_DOCS "Standard C++ minimal d'une bibliothèque wsl")

wsl_add_library(array2d      STANDARD 23)
wsl_add_library(dlist        STANDARD 23)
wsl_add_library(graded       STANDARD 23)
wsl_add_library(graded_io    STANDARD 23 DEPENDS graded)
wsl_add_library(elo          STANDARD 23)
wsl_add_library(progress     STANDARD 23)
wsl_add_library(special      STANDARD 23)
wsl_add_library(permutation  STANDARD 23)
wsl_add_library(random       STANDARD 23 DEPENDS special permutation)

# Bibliothèques dont on construit les tests et exemples (toutes, sauf en CI).
set(WSL_LIBRARIES_TESTED ALL)
set(WSL_LIBRARIES_TESTED_EXAMPLES "")
if(WSL_LIBRARIES_AFFECTED_BY)
    include(cmake/affected.cmake)
    wsl_affected_libraries(${WSL_LIBRARIES_AFFECTED_BY} WSL_LIBRARIES_TESTED)
    message(STATUS "Bibliothèques testées : ${WSL_LIBRARIES_TESTED}")
endif()

# Vrai si les tests de <bibliothèque> sont à construire.
function(wsl_is_tested library variable)
    if(WSL_LIBRARIES_TESTED STREQUAL "ALL" OR library IN_LIST WSL_LIBRARIES_TESTED)
        set(${variable} TRUE PARENT_SCOPE)
    else()
        set(${variable} FALSE PARENT_SCOPE)
    endif()
endfunction()

# Tout d'un coup : wsl::libraries.
add_library(wsl_libraries INTERFACE)
add_library(wsl::libraries ALIAS wsl_libraries)
set_target_properties(wsl_libraries PROPERTIES EXPORT_NAME libraries)
target_link_libraries(wsl_libraries INTERFACE ${WSL_LIBRARIES_ALL})

# ---------------------------------------------------------------------------
# Options de compilation des tests et exemples (jamais transmises aux consommateurs)
# ---------------------------------------------------------------------------
add_library(wsl_libraries_options INTERFACE)
target_compile_options(wsl_libraries_options INTERFACE -Wall -Wextra -Wpedantic)
if(WSL_LIBRARIES_WERROR)
    target_compile_options(wsl_libraries_options INTERFACE -Werror)
endif()
if(WSL_LIBRARIES_SANITIZERS)
    target_compile_options(wsl_libraries_options INTERFACE
        -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(wsl_libraries_options INTERFACE -fsanitize=address,undefined)
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
#
#     wsl_add_test(<nom> [LIBRARY <bibliothèque>])
#
# tests/test_<nom>.cpp, lié à wsl::<bibliothèque> (par défaut wsl::<nom>), et
# compilé une fois par standard supporté par cette bibliothèque :
# test_<nom> en C++23 (celui que lancent l'éditeur et le débogueur), et
# test_<nom>_cxx<standard> pour chaque standard inférieur supporté. Un test
# réussit s'il renvoie 0.

if(WSL_LIBRARIES_BUILD_TESTS)
    enable_testing()

    function(wsl_add_test name)
        cmake_parse_arguments(PARSE_ARGV 1 ARG "" "LIBRARY" "")
        set(library ${name})
        if(ARG_LIBRARY)
            set(library ${ARG_LIBRARY})
        endif()
        wsl_is_tested(${library} tested)
        if(NOT tested)
            return()
        endif()
        get_target_property(minimum wsl_${library} WSL_STANDARD)
        foreach(standard IN LISTS WSL_LIBRARIES_STANDARDS)
            if(standard LESS minimum)
                continue()
            endif()
            if(standard EQUAL 23)
                set(test test_${name})
            else()
                set(test test_${name}_cxx${standard})
            endif()
            add_executable(${test} tests/test_${name}.cpp)
            target_link_libraries(${test} PRIVATE wsl::${library} wsl_libraries_options)
            set_target_properties(${test} PROPERTIES
                CXX_STANDARD ${standard}
                CXX_STANDARD_REQUIRED ON
                RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
            # Les assert restent actifs, même en Release (qui définit NDEBUG).
            target_compile_options(${test} PRIVATE -UNDEBUG)
            add_test(NAME ${test} COMMAND ${test})
        endforeach()
    endfunction()

    wsl_add_test(array2d)
    wsl_add_test(dlist)
    wsl_add_test(graded)
    # Le test compare exprès jet<int> et jet<unsigned>, pour vérifier que le
    # piège signé/non signé des scalaires est reproduit tel quel : l'avertissement
    # qu'en reçoit un utilisateur est voulu, mais n'a pas sa place ici.
    # (Par une cible liée après wsl_libraries_options : une option propre à la
    # cible passerait avant -Wextra, qui la réactiverait.)
    add_library(wsl_libraries_no_sign_compare INTERFACE)
    target_compile_options(wsl_libraries_no_sign_compare INTERFACE -Wno-sign-compare)
    foreach(test test_graded test_graded_cxx20)
        if(NOT TARGET ${test})
            continue()
        endif()
        target_link_libraries(${test} PRIVATE wsl_libraries_no_sign_compare)
    endforeach()
    wsl_add_test(graded_io)
    wsl_add_test(elo)
    wsl_add_test(progress)
    wsl_add_test(permutation)
    # Valeurs de référence de wsl::special (mpmath) : engendrées dans le dossier de
    # build, pas versionnées. Il faut uv, ou un python3 muni de mpmath ; sinon, seul
    # le test de précision est sauté.
    find_program(WSL_UV uv HINTS $ENV{HOME}/.local/bin)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    set(reference_generator "")
    if(WSL_UV)
        set(reference_generator ${WSL_UV} run --no-project --quiet --with mpmath python)
    elseif(Python3_Interpreter_FOUND)
        execute_process(COMMAND ${Python3_EXECUTABLE} -c "import mpmath"
                        RESULT_VARIABLE no_mpmath OUTPUT_QUIET ERROR_QUIET)
        if(no_mpmath EQUAL 0)
            set(reference_generator ${Python3_EXECUTABLE})
        endif()
    endif()

    wsl_is_tested(special special_tested)
    wsl_is_tested(random random_tested)
    if(reference_generator AND special_tested)
        set(references ${PROJECT_BINARY_DIR}/generated/data/special_elementary_reference.hpp)
        add_custom_command(
            OUTPUT ${references}
            COMMAND ${reference_generator} ${PROJECT_SOURCE_DIR}/tools/special/generate.py --references ${references}
            DEPENDS ${PROJECT_SOURCE_DIR}/tools/special/generate.py
            COMMENT "Valeurs de référence de wsl::special (mpmath)")
        add_custom_target(special_references DEPENDS ${references})
        wsl_add_test(special_elementary LIBRARY special)
        foreach(test test_special_elementary test_special_elementary_cxx20)
            if(NOT TARGET ${test})
                continue()
            endif()
            add_dependencies(${test} special_references)
            target_include_directories(${test} PRIVATE ${PROJECT_BINARY_DIR}/generated)
        endforeach()
    endif()
    if(reference_generator AND random_tested)
        # Lois de wsl::random : même principe, une famille par module de tools/random/references/
        # (seules les familles modifiées sont recalculées) ; l'en-tête commun les inclut toutes.
        set(law_families base derived circular complement inverse_gaussian skew_normal stable discrete skellam zipf rice)
        set(law_sources ${PROJECT_SOURCE_DIR}/tools/random/references)
        set(law_outputs "")
        set(law_includes "")
        foreach(family IN LISTS law_families)
            set(output ${PROJECT_BINARY_DIR}/generated/data/random_law_reference/${family}.hpp)
            set(extra "")
            if(family STREQUAL "base")
                set(extra ${PROJECT_SOURCE_DIR}/include/wsl/special/detail/extended.hpp)
            elseif(family STREQUAL "derived")
                set(extra ${law_sources}/base.py)
            endif()
            add_custom_command(
                OUTPUT ${output}
                COMMAND ${reference_generator} ${PROJECT_SOURCE_DIR}/tools/random/law_references.py ${family} ${output}
                DEPENDS ${PROJECT_SOURCE_DIR}/tools/random/law_references.py ${law_sources}/${family}.py ${law_sources}/common.py
                        ${PROJECT_SOURCE_DIR}/tools/special/generate.py ${extra}
                COMMENT "Valeurs de référence des lois de wsl::random : ${family} (mpmath)")
            list(APPEND law_outputs ${output})
            string(APPEND law_includes "#include \"data/random_law_reference/${family}.hpp\"\n")
        endforeach()
        file(CONFIGURE OUTPUT ${PROJECT_BINARY_DIR}/generated/data/random_law_reference.hpp CONTENT [[
// ENGENDRÉ par CMake — ne pas modifier à la main.
#pragma once

// Valeurs de référence des lois de wsl::random (mpmath), une famille par en-tête.

namespace law_reference_data
{
    template <class T>
    struct point
    {
        double a, b, x;
        T rounded;
        double offset;
    };

    struct quantile_point
    {
        double p, u;
        long long k;
    };

    struct quantile_point2
    {
        double a, b, u;
        long long k;
    };
}

@law_includes@]] @ONLY)
        add_custom_target(random_law_references DEPENDS ${law_outputs})
        wsl_add_test(random_precision LIBRARY random)
        foreach(test test_random_precision test_random_precision_cxx20)
            if(NOT TARGET ${test})
                continue()
            endif()
            add_dependencies(${test} random_law_references)
            target_include_directories(${test} PRIVATE ${PROJECT_BINARY_DIR}/generated)
        endforeach()
    endif()
    if(NOT reference_generator)
        message(WARNING "Ni uv ni python3 avec mpmath : tests de précision de wsl::special et wsl::random sautés")
    endif()
    wsl_add_test(special_bitexact LIBRARY special)
    wsl_add_test(random_engines LIBRARY random)
    wsl_add_test(random_uniform LIBRARY random)
    wsl_add_test(random_laws LIBRARY random)
    wsl_add_test(random_algebra LIBRARY random)
    wsl_add_test(random_multivariate LIBRARY random)
    wsl_add_test(random_combinatorics LIBRARY random)
    wsl_add_test(random_processes LIBRARY random)
    wsl_add_test(random_continuous LIBRARY random)
    wsl_add_test(random_variables LIBRARY random)
    # Aucun produit de wsl::special ne doit pouvoir être fusionné en fma.
    if(special_tested)
        add_test(NAME special_no_contraction
                 COMMAND ${CMAKE_COMMAND} -DCXX=${CMAKE_CXX_COMPILER} -DINCLUDE=${PROJECT_SOURCE_DIR}/include
                         -DSOURCES=${PROJECT_SOURCE_DIR}/tests/special_contraction
                         -DPROCESSOR=${CMAKE_SYSTEM_PROCESSOR}
                         -P ${PROJECT_SOURCE_DIR}/cmake/check_no_contraction.cmake)
    endif()
endif()

# ---------------------------------------------------------------------------
# Exemples : examples/<nom>.cpp, un exécutable chacun
# ---------------------------------------------------------------------------
if(WSL_LIBRARIES_BUILD_EXAMPLES)
    file(GLOB example_sources CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/examples/*.cpp)
    foreach(source IN LISTS example_sources)
        get_filename_component(example ${source} NAME_WE)
        # En CI, un exemple n'est construit que s'il a changé ou s'il inclut une
        # bibliothèque testée.
        if(NOT WSL_LIBRARIES_TESTED STREQUAL "ALL" AND NOT example IN_LIST WSL_LIBRARIES_TESTED_EXAMPLES)
            file(STRINGS ${source} includes REGEX "^#include <wsl/[a-z0-9_]+")
            set(wanted FALSE)
            foreach(line IN LISTS includes)
                string(REGEX REPLACE "^#include <wsl/([a-z0-9_]+).*" "\\1" library "${line}")
                if(library IN_LIST WSL_LIBRARIES_TESTED)
                    set(wanted TRUE)
                endif()
            endforeach()
            if(NOT wanted)
                continue()
            endif()
        endif()
        add_executable(${example} ${source})
        target_link_libraries(${example} PRIVATE wsl::libraries wsl_libraries_options)
        set_target_properties(${example} PROPERTIES
            CXX_STANDARD 23
            RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/examples)
    endforeach()
endif()

# ---------------------------------------------------------------------------
# Installation : find_package(wsl_libraries) chez qui nous utilise
# ---------------------------------------------------------------------------
#
#     cmake --install build/release --prefix ~/.local
#
# En-têtes seuls : rien ne dépend de l'architecture, d'où share/ plutôt que lib/.

include(CMakePackageConfigHelpers)

set(WSL_LIBRARIES_CMAKE_DIR ${CMAKE_INSTALL_DATADIR}/cmake/wsl_libraries)

install(TARGETS ${WSL_LIBRARIES_ALL} wsl_libraries
        EXPORT wsl_librariesTargets
        FILE_SET HEADERS)

install(EXPORT wsl_librariesTargets
        NAMESPACE wsl::
        DESTINATION ${WSL_LIBRARIES_CMAKE_DIR})

configure_package_config_file(
    cmake/wsl_librariesConfig.cmake.in
    ${PROJECT_BINARY_DIR}/wsl_librariesConfig.cmake
    INSTALL_DESTINATION ${WSL_LIBRARIES_CMAKE_DIR})

# Avant 1.0, une version mineure peut casser l'interface : seule la même
# version mineure est déclarée compatible.
write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/wsl_librariesConfigVersion.cmake
    COMPATIBILITY SameMinorVersion
    ARCH_INDEPENDENT)

install(FILES
    ${PROJECT_BINARY_DIR}/wsl_librariesConfig.cmake
    ${PROJECT_BINARY_DIR}/wsl_librariesConfigVersion.cmake
    DESTINATION ${WSL_LIBRARIES_CMAKE_DIR})
