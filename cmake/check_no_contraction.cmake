# Vérifie que wsl::special ne laisse fusionner aucun produit en fma.
#
#     cmake -DCXX=<compilateur> -DINCLUDE=<include/> -DSOURCES=<dossier> -DPROCESSOR=<arch> -P check_no_contraction.cmake
#
# Compile tests/special_contraction/functions.cpp en assembleur, fma disponible et
# contraction permise partout, et y cherche une instruction fma. Contre-épreuve :
# control.cpp, compilé pareil, doit en contenir une.
#
# Limite : les compilateurs ne fusionnent un produit que s'il n'a qu'un usage. Un
# produit non protégé mais utilisé plusieurs fois échappe donc à ce test (vérifié) ;
# il ne change aucun résultat aujourd'hui, mais un compilateur futur pourrait le
# fusionner : les empreintes de test_special_bitexact le verraient alors.

if(PROCESSOR MATCHES "x86_64|AMD64|amd64")
    set(target_flags -march=x86-64-v3)
elseif(PROCESSOR MATCHES "arm64|aarch64")
    set(target_flags "")
else()
    message(STATUS "Architecture ${PROCESSOR} : vérification non prévue, ignorée")
    return()
endif()

set(fma_regex "[ \t](v?fn?m(add|sub)[0-9a-z]*)[ \t]")

function(assembly source out)
    execute_process(
        COMMAND ${CXX} -std=c++20 -O2 -ffp-contract=fast ${target_flags} -I${INCLUDE} -S -o - ${source}
        OUTPUT_VARIABLE asm ERROR_VARIABLE errors RESULT_VARIABLE status)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR "Compilation de ${source} impossible :\n${errors}")
    endif()
    set(${out} "${asm}" PARENT_SCOPE)
endfunction()

assembly(${SOURCES}/control.cpp control)
if(NOT control MATCHES "${fma_regex}")
    message(FATAL_ERROR "Contre-épreuve : aucune fma dans control.cpp, la vérification ne prouverait rien")
endif()
set(control_instruction ${CMAKE_MATCH_1})

assembly(${SOURCES}/functions.cpp functions)
string(REGEX MATCHALL "${fma_regex}" found "${functions}")
if(found)
    message(FATAL_ERROR "fma dans le code de wsl::special : ${found}\n"
                        "Un produit entrant dans une addition doit passer par detail::fmul.")
endif()
message(STATUS "Aucune fma dans wsl::special (contre-épreuve : ${control_instruction})")
