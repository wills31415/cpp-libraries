# Tests qui doivent échouer à la compilation, avec le message attendu.
#
#     cmake -DCXX=<compilateur> -DINCLUDE=<include> -DSOURCE=<fichier.cpp> -P check_compile_fails.cmake
#
# Le fichier contient des cas `#if CASE == <n>` et, pour chacun, une ligne
#     // ERREUR <n> : <texte attendu dans la sortie du compilateur>
# Sans CASE, le fichier doit compiler (sinon un cas échouerait pour une autre raison).

file(STRINGS ${SOURCE} lines ENCODING UTF-8 REGEX "^// ERREUR [0-9]+ : ")
if(NOT lines)
    message(FATAL_ERROR "${SOURCE} : aucun cas « // ERREUR <n> : <texte> »")
endif()

set(common -std=c++23 -fsyntax-only -I${INCLUDE} ${SOURCE})

execute_process(COMMAND ${CXX} ${common} RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE output)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "${SOURCE} ne compile pas sans CASE :\n${output}")
endif()

set(failures 0)
foreach(line IN LISTS lines)
    # Le texte attendu est ce qui suit « : » (fichier lu en UTF-8 : sans ENCODING,
    # file(STRINGS) coupe les lignes aux caractères accentués).
    string(REGEX MATCH "^// ERREUR ([0-9]+) : " prefix "${line}")
    set(case ${CMAKE_MATCH_1})
    string(LENGTH "${prefix}" length)
    string(SUBSTRING "${line}" ${length} -1 expected)
    string(STRIP "${expected}" expected)
    if(expected STREQUAL "")
        message(FATAL_ERROR "cas ${case} : texte attendu vide")
    endif()
    execute_process(COMMAND ${CXX} ${common} -DCASE=${case}
                    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE output)
    if(status EQUAL 0)
        message(SEND_ERROR "cas ${case} : compile, alors qu'il devait être refusé (${expected})")
        math(EXPR failures "${failures} + 1")
    else()
        string(FIND "${output}" "${expected}" position)
        if(position EQUAL -1)
            message(SEND_ERROR "cas ${case} : refusé, mais sans le message attendu « ${expected} » :\n${output}")
            math(EXPR failures "${failures} + 1")
        else()
            message(STATUS "cas ${case} : refusé (${expected})")
        endif()
    endif()
endforeach()

if(failures GREATER 0)
    message(FATAL_ERROR "${failures} cas incorrects dans ${SOURCE}")
endif()
