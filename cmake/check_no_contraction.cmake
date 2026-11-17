# Vérifie que wsl::special ne laisse fusionner aucun produit en fma.
#
#     cmake -DCXX=<compilateur> -DINCLUDE=<include/> -DSOURCES=<dossier> -DPROCESSOR=<arch> -P check_no_contraction.cmake
#
# Compile tests/special_contraction/functions.cpp en assembleur, fma disponible et
# contraction permise partout, et y cherche une instruction fma. Contre-épreuve :
# control.cpp, compilé pareil, doit en contenir une. Les fma explicites des
# produits exacts (two_product, quand la cible en a) sont écartées par
# WSL_SPECIAL_NO_FMA : seules restent visibles les contractions involontaires.
#
# Limite : les compilateurs ne fusionnent un produit que s'il n'a qu'un usage. Un
# produit non protégé mais utilisé plusieurs fois échappe donc à ce test (vérifié) ;
# il ne change aucun résultat aujourd'hui, mais un compilateur futur pourrait le
# fusionner : les empreintes de test_special_bitexact le verraient alors.

# Les deux architectures ne fusionnent pas les mêmes produits (ARM64 fusionne aussi
# des produits à usages multiples, que x86-64 laisse) : avec Apple Clang, l'autre
# architecture est vérifiée en plus, par compilation croisée (assembleur seulement).

if(PROCESSOR MATCHES "x86_64|AMD64|amd64")
    set(configurations "native:-march=x86-64-v3")
    set(other "arm64-apple-macos:")
elseif(PROCESSOR MATCHES "arm64|aarch64")
    set(configurations "native:")
    set(other "x86_64-apple-macos:-march=x86-64-v3")
else()
    message(STATUS "Architecture ${PROCESSOR} : vérification non prévue, ignorée")
    return()
endif()
execute_process(COMMAND ${CXX} --version OUTPUT_VARIABLE version)
if(version MATCHES "Apple clang")
    list(APPEND configurations ${other})
endif()

set(fma_regex "[ \t](v?fn?m(add|sub)[0-9a-z]*)[ \t]")

function(assembly source target flags out)
    set(target_option "")
    if(NOT target STREQUAL "native")
        set(target_option -target ${target})
    endif()
    execute_process(
        COMMAND ${CXX} -std=c++20 -O2 -ffp-contract=fast -DWSL_SPECIAL_NO_FMA ${target_option} ${flags} -I${INCLUDE} -S -o - ${source}
        OUTPUT_VARIABLE asm ERROR_VARIABLE errors RESULT_VARIABLE status)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR "Compilation de ${source} (${target}) impossible :\n${errors}")
    endif()
    set(${out} "${asm}" PARENT_SCOPE)
endfunction()

foreach(configuration IN LISTS configurations)
    string(REPLACE ":" ";" parts "${configuration}:")
    list(GET parts 0 target)
    list(GET parts 1 flags)
    separate_arguments(flags)

    assembly(${SOURCES}/control.cpp "${target}" "${flags}" control)
    if(NOT control MATCHES "${fma_regex}")
        message(FATAL_ERROR "Contre-épreuve (${target}) : aucune fma dans control.cpp, la vérification ne prouverait rien")
    endif()
    set(control_instruction ${CMAKE_MATCH_1})

    assembly(${SOURCES}/functions.cpp "${target}" "${flags}" functions)
    string(REGEX MATCHALL "${fma_regex}" found "${functions}")
    if(found)
        message(FATAL_ERROR "fma dans le code de wsl::special (${target}) : ${found}\n"
                            "Un produit entrant dans une addition doit passer par detail::fmul.")
    endif()
    message(STATUS "${target} : aucune fma dans wsl::special (contre-épreuve : ${control_instruction})")
endforeach()
