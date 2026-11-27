cmake_minimum_required(VERSION 3.28)

project(wsl
        VERSION 0.1.0
        DESCRIPTION "Bibliothèque scientifique C++23 en en-têtes seuls (nom provisoire)"
        LANGUAGES CXX)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Tests et exemples : construits seulement quand ce dépôt est le projet principal.
# Un projet qui nous intègre par FetchContent ou add_subdirectory ne les voit pas.
option(WSL_BUILD_TESTS "Construire les tests" ${PROJECT_IS_TOP_LEVEL})
option(WSL_BUILD_EXAMPLES "Construire les exemples" ${PROJECT_IS_TOP_LEVEL})
option(WSL_SANITIZERS "Compiler tests et exemples avec ASan et UBSan" OFF)
option(WSL_WERROR "Traiter les avertissements des tests et exemples comme des erreurs (CI)" OFF)
# CI : fichier listant les chemins modifiés (un par ligne, « ALL » pour tout) ; seuls
# les tests et exemples des bibliothèques touchées, et de celles qui en dépendent, sont
# construits (cmake/affected.cmake). Vide : tout.
set(WSL_AFFECTED_BY "" CACHE FILEPATH "Fichier des chemins modifiés (CI)")

include(GNUInstallDirs)

# Standards sur lesquels les tests sont compilés. Une bibliothèque est testée
# sur chacun de ceux qu'elle prétend supporter (de son minimum à C++23). Le dépôt est
# en C++23 partout (design/interfaces.md, décision 7) ; ajouter 20 ici si une
# bibliothèque redescend un jour à C++20, sur décision explicite.
set(WSL_STANDARDS 23)

# ---------------------------------------------------------------------------
# Les bibliothèques
# ---------------------------------------------------------------------------
#
#     wsl_add_library(<nom> STANDARD <20|23> [DEPENDS <nom>...])
#
# <nom> est l'espace de noms sous wsl:: (design/interfaces.md, décision 31), avec ses
# « :: » : random, random::engines, physics::continuum_mechanics::fluids… La cible est
# wsl::<nom> ; ses en-têtes sont include/wsl/<chemin>.hpp et include/wsl/<chemin>/ (le
# chemin est le nom, « :: » devenant « / »), sauf ceux d'une bibliothèque plus profonde,
# qu'elle inclut (une cible parent inclut ses enfants). STANDARD est le minimum exigé,
# transmis aux consommateurs ; DEPENDS nomme d'autres bibliothèques du dépôt.
# wsl_finalize_libraries(), après la dernière déclaration, attribue les en-têtes.

set(WSL_ALL "")

# En-têtes sans extension : relais engendrés à la configuration (cmake/relays.cmake).
include(cmake/relays.cmake)
wsl_generate_relays(${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/relays WSL_RELAY_ROOTS)

# Identifiant CMake d'une bibliothèque : « :: » est réservé aux cibles importées et alias.
function(wsl_library_id name variable)
    string(REPLACE "::" "__" id "${name}")
    set(${variable} ${id} PARENT_SCOPE)
endfunction()

function(wsl_add_library name)
    cmake_parse_arguments(PARSE_ARGV 1 ARG "" "STANDARD" "DEPENDS")
    if(NOT ARG_STANDARD)
        message(FATAL_ERROR "wsl_add_library(${name}) : STANDARD manquant")
    endif()

    wsl_library_id(${name} id)
    string(REPLACE "::" "/" path "${name}")
    set(target wsl_${id})
    add_library(${target} INTERFACE)
    add_library(wsl::${name} ALIAS ${target})
    set_target_properties(${target} PROPERTIES
        EXPORT_NAME ${name}
        WSL_STANDARD ${ARG_STANDARD})
    target_compile_features(${target} INTERFACE cxx_std_${ARG_STANDARD})
    foreach(root IN LISTS WSL_RELAY_ROOTS)
        target_include_directories(${target} INTERFACE
            $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/relays/${root}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/wsl-relays/${root}>)
    endforeach()
    foreach(dependency IN LISTS ARG_DEPENDS)
        wsl_library_id(${dependency} dependency_id)
        target_link_libraries(${target} INTERFACE wsl_${dependency_id})
    endforeach()

    set_property(GLOBAL APPEND PROPERTY WSL_LIBRARY_NAMES ${name})
    set_property(GLOBAL PROPERTY WSL_PATH_${id} ${path})
    set_property(GLOBAL PROPERTY WSL_DEPENDS_${id} ${ARG_DEPENDS})
    set(WSL_ALL ${WSL_ALL} ${target} PARENT_SCOPE)
endfunction()

# Bibliothèque la plus profonde dont <chemin> (relatif à la racine d'un arbre : wsl/… pour
# include/, sans wsl/ pour tests/) relève ; vide si aucune.
function(wsl_owner relative variable)
    get_property(names GLOBAL PROPERTY WSL_LIBRARY_NAMES)
    set(best "")
    set(best_length -1)
    foreach(name IN LISTS names)
        wsl_library_id(${name} id)
        get_property(path GLOBAL PROPERTY WSL_PATH_${id})
        string(LENGTH "${path}" length)
        if(relative MATCHES "^${path}(\\.[a-z]+|/|$)" AND length GREATER best_length)
            set(best ${name})
            set(best_length ${length})
        endif()
    endforeach()
    set(${variable} "${best}" PARENT_SCOPE)
endfunction()

function(wsl_finalize_libraries)
    get_property(names GLOBAL PROPERTY WSL_LIBRARY_NAMES)
    file(GLOB_RECURSE headers CONFIGURE_DEPENDS RELATIVE ${PROJECT_SOURCE_DIR}/include
         ${PROJECT_SOURCE_DIR}/include/wsl/*.hpp)
    foreach(header IN LISTS headers)
        string(REGEX REPLACE "^wsl/" "" relative ${header})
        wsl_owner(${relative} owner)
        if(owner STREQUAL "")
            message(FATAL_ERROR "include/${header} n'appartient à aucune bibliothèque (wsl_add_library)")
        endif()
        wsl_library_id(${owner} id)
        set_property(GLOBAL APPEND PROPERTY WSL_HEADERS_${id} ${header})
    endforeach()
    set(data "")
    foreach(name IN LISTS names)
        wsl_library_id(${name} id)
        get_property(own GLOBAL PROPERTY WSL_HEADERS_${id})
        list(TRANSFORM own PREPEND ${PROJECT_SOURCE_DIR}/include/)
        target_sources(wsl_${id} INTERFACE
            FILE_SET HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include FILES ${own})
        # Une cible parent inclut ses enfants.
        get_property(path GLOBAL PROPERTY WSL_PATH_${id})
        foreach(other IN LISTS names)
            wsl_library_id(${other} other_id)
            get_property(other_path GLOBAL PROPERTY WSL_PATH_${other_id})
            if(other_path MATCHES "^${path}/")
                target_link_libraries(wsl_${id} INTERFACE wsl_${other_id})
            endif()
        endforeach()
        get_property(dependencies GLOBAL PROPERTY WSL_DEPENDS_${id})
        string(REPLACE ";" "," dependencies "${dependencies}")
        get_property(own GLOBAL PROPERTY WSL_HEADERS_${id})
        string(REPLACE ";" "," own "${own}")
        string(APPEND data "${name}|${path}|${dependencies}|${own}\n")
    endforeach()
    # Pour le contrôle des inclusions (cmake/check_includes.cmake).
    file(CONFIGURE OUTPUT ${PROJECT_BINARY_DIR}/libraries.txt CONTENT "${data}")
endfunction()

define_property(TARGET PROPERTY WSL_STANDARD
    BRIEF_DOCS "Standard C++ minimal d'une bibliothèque wsl")

# Anciennes bibliothèques (array2d, graded, graded_io, elo, progress) : dans legacy/,
# figées, ni construites ni installées, en attendant leur refonte (design/interfaces.md,
# décision 36).
wsl_add_library(special      STANDARD 23)
wsl_add_library(permutation  STANDARD 23)
wsl_add_library(random       STANDARD 23 DEPENDS special permutation)
wsl_add_library(mdarray      STANDARD 23 DEPENDS special)
wsl_finalize_libraries()
# Mode de calcul flottant des instructions non marquées de wsl::mdarray (design/mdarray.md,
# I9) : rapide par défaut ; reproductible au bit près partout si ON. Transmis à tout le
# projet consommateur (une même valeur partout : règle de définition unique).
option(WSL_MDARRAY_REPRODUCIBLE "wsl::mdarray : mode reproductible par défaut" OFF)
if(WSL_MDARRAY_REPRODUCIBLE)
    target_compile_definitions(wsl_mdarray INTERFACE WSL_MDARRAY_REPRODUCIBLE=1)
endif()

# Bibliothèques dont on construit les tests et exemples (toutes, sauf en CI).
set(WSL_TESTED ALL)
set(WSL_TESTED_EXAMPLES "")
if(WSL_AFFECTED_BY)
    include(cmake/affected.cmake)
    wsl_affected_libraries(${WSL_AFFECTED_BY} WSL_TESTED)
    message(STATUS "Bibliothèques testées : ${WSL_TESTED}")
endif()

# Vrai si les tests de <bibliothèque> sont à construire.
function(wsl_is_tested library variable)
    if(WSL_TESTED STREQUAL "ALL" OR library IN_LIST WSL_TESTED)
        set(${variable} TRUE PARENT_SCOPE)
    else()
        set(${variable} FALSE PARENT_SCOPE)
    endif()
endfunction()

# Tout d'un coup : wsl::wsl (design/interfaces.md, décision 31).
add_library(wsl_all INTERFACE)
add_library(wsl::wsl ALIAS wsl_all)
set_target_properties(wsl_all PROPERTIES EXPORT_NAME wsl)
target_link_libraries(wsl_all INTERFACE ${WSL_ALL})

# ---------------------------------------------------------------------------
# Options de compilation des tests et exemples (jamais transmises aux consommateurs)
# ---------------------------------------------------------------------------
add_library(wsl_options INTERFACE)
target_compile_options(wsl_options INTERFACE -Wall -Wextra -Wpedantic)
if(WSL_WERROR)
    target_compile_options(wsl_options INTERFACE -Werror)
endif()
if(WSL_SANITIZERS)
    target_compile_options(wsl_options INTERFACE
        -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(wsl_options INTERFACE -fsanitize=address,undefined)
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

if(WSL_BUILD_TESTS)
    enable_testing()

    function(wsl_add_test name)
        cmake_parse_arguments(PARSE_ARGV 1 ARG "" "LIBRARY" "LABELS")
        if(NOT ARG_LABELS)
            set(ARG_LABELS properties)
        endif()
        set(library ${name})
        if(ARG_LIBRARY)
            set(library ${ARG_LIBRARY})
        endif()
        wsl_is_tested(${library} tested)
        if(NOT tested)
            return()
        endif()
        get_target_property(minimum wsl_${library} WSL_STANDARD)
        foreach(standard IN LISTS WSL_STANDARDS)
            if(standard LESS minimum)
                continue()
            endif()
            if(standard EQUAL 23)
                set(test test_${name})
            else()
                set(test test_${name}_cxx${standard})
            endif()
            add_executable(${test} tests/test_${name}.cpp)
            target_link_libraries(${test} PRIVATE wsl::${library} wsl_options)
            set_target_properties(${test} PROPERTIES
                CXX_STANDARD ${standard}
                CXX_STANDARD_REQUIRED ON
                RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
            # Les assert restent actifs, même en Release (qui définit NDEBUG).
            target_compile_options(${test} PRIVATE -UNDEBUG)
            add_test(NAME ${test} COMMAND ${test})
            set_tests_properties(${test} PROPERTIES LABELS "${ARG_LABELS}")
        endforeach()
    endfunction()

    # Tests rangés comme les espaces de noms (design/interfaces.md, décision 32) :
    # tests/<chemin>.cpp, nommé d'après son chemin (tests/mdarray/axis.cpp →
    # « mdarray.axis »), lié à la bibliothèque la plus profonde dont il relève. Un dossier
    # compile_fail/ contient des cas qui doivent échouer à la compilation avec le message
    # attendu (cmake/check_compile_fails.cmake). Étiquettes CTest : la nature de la batterie
    # et son coût, déclarés en tête du fichier par une ligne
    #     // Étiquettes : precision long
    # (par défaut « properties » ; « compile_fail » pour un cas de refus). Les tests de
    # special et random (tests/test_*.cpp) suivront avec leur découpage.
    file(GLOB_RECURSE tree_tests CONFIGURE_DEPENDS RELATIVE ${PROJECT_SOURCE_DIR}/tests
         ${PROJECT_SOURCE_DIR}/tests/*.cpp)
    foreach(relative IN LISTS tree_tests)
        if(relative MATCHES "^(test_|package/|special_contraction/)")
            continue()
        endif()
        string(REGEX REPLACE "\\.cpp$" "" stem ${relative})
        string(REPLACE "/" "." test ${stem})
        wsl_owner(${stem} library)
        if(library STREQUAL "")
            message(FATAL_ERROR "tests/${relative} ne relève d'aucune bibliothèque")
        endif()
        wsl_is_tested(${library} tested)
        if(NOT tested)
            continue()
        endif()
        set(source ${PROJECT_SOURCE_DIR}/tests/${relative})
        file(STRINGS ${source} declared LIMIT_COUNT 1 ENCODING UTF-8 REGEX "^// Étiquettes : ")
        string(REGEX REPLACE "^// Étiquettes : " "" labels "${declared}")
        separate_arguments(labels)
        if(stem MATCHES "(^|/)compile_fail/")
            add_test(NAME ${test}
                     COMMAND ${CMAKE_COMMAND} -DCXX=${CMAKE_CXX_COMPILER} -DINCLUDE=${PROJECT_SOURCE_DIR}/include
                             -DSOURCE=${source} -P ${PROJECT_SOURCE_DIR}/cmake/check_compile_fails.cmake)
            list(APPEND labels compile_fail)
        else()
            add_executable(${test} ${source})
            wsl_library_id(${library} id)
            target_link_libraries(${test} PRIVATE wsl_${id} wsl_options)
            set_target_properties(${test} PROPERTIES
                CXX_STANDARD 23
                CXX_STANDARD_REQUIRED ON
                RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
            target_compile_options(${test} PRIVATE -UNDEBUG)
            add_test(NAME ${test} COMMAND ${test})
            if(NOT labels)
                set(labels properties)
            endif()
        endif()
        set_tests_properties(${test} PROPERTIES LABELS "${labels}")
    endforeach()

    # Contrôle des inclusions (décision 31) : chaque #include <wsl/…> d'un en-tête relève
    # de sa bibliothèque, d'un de ses enfants ou d'une de ses dépendances.
    add_test(NAME infrastructure.includes
             COMMAND ${CMAKE_COMMAND} -DDATA=${PROJECT_BINARY_DIR}/libraries.txt
                     -DINCLUDE=${PROJECT_SOURCE_DIR}/include -P ${PROJECT_SOURCE_DIR}/cmake/check_includes.cmake)
    set_tests_properties(infrastructure.includes PROPERTIES LABELS properties)

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
        wsl_add_test(special_elementary LIBRARY special LABELS precision long)
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
        wsl_add_test(random_precision LIBRARY random LABELS precision long)
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
    wsl_add_test(special_bitexact LIBRARY special LABELS fingerprints)
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
if(WSL_BUILD_EXAMPLES)
    file(GLOB example_sources CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/examples/*.cpp)
    foreach(source IN LISTS example_sources)
        get_filename_component(example ${source} NAME_WE)
        # En CI, un exemple n'est construit que s'il a changé ou s'il inclut une
        # bibliothèque testée.
        if(NOT WSL_TESTED STREQUAL "ALL" AND NOT example IN_LIST WSL_TESTED_EXAMPLES)
            file(STRINGS ${source} includes REGEX "^#include <wsl/[a-z0-9_]+")
            set(wanted FALSE)
            foreach(line IN LISTS includes)
                string(REGEX REPLACE "^#include <wsl/([a-z0-9_]+).*" "\\1" library "${line}")
                if(library IN_LIST WSL_TESTED)
                    set(wanted TRUE)
                endif()
            endforeach()
            if(NOT wanted)
                continue()
            endif()
        endif()
        add_executable(${example} ${source})
        target_link_libraries(${example} PRIVATE wsl::wsl wsl_options)
        set_target_properties(${example} PROPERTIES
            CXX_STANDARD 23
            RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/examples)
    endforeach()
endif()

# ---------------------------------------------------------------------------
# Installation : find_package(wsl) chez qui nous utilise
# ---------------------------------------------------------------------------
#
#     cmake --install build/release --prefix ~/.local
#
# En-têtes seuls : rien ne dépend de l'architecture, d'où share/ plutôt que lib/.

include(CMakePackageConfigHelpers)

set(WSL_CMAKE_DIR ${CMAKE_INSTALL_DATADIR}/cmake/wsl)

install(TARGETS ${WSL_ALL} wsl_all
        EXPORT wslTargets
        FILE_SET HEADERS)

install(DIRECTORY ${PROJECT_BINARY_DIR}/relays/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wsl-relays)

install(EXPORT wslTargets
        NAMESPACE wsl::
        DESTINATION ${WSL_CMAKE_DIR})

configure_package_config_file(
    cmake/wslConfig.cmake.in
    ${PROJECT_BINARY_DIR}/wslConfig.cmake
    INSTALL_DESTINATION ${WSL_CMAKE_DIR})

# Avant 1.0, une version mineure peut casser l'interface : seule la même
# version mineure est déclarée compatible.
write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/wslConfigVersion.cmake
    COMPATIBILITY SameMinorVersion
    ARCH_INDEPENDENT)

install(FILES
    ${PROJECT_BINARY_DIR}/wslConfig.cmake
    ${PROJECT_BINARY_DIR}/wslConfigVersion.cmake
    DESTINATION ${WSL_CMAKE_DIR})
