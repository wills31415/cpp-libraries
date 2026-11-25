# Bibliothèques touchées par une liste de fichiers modifiés (CI : ne compiler et ne
# tester que ce qu'un push peut avoir changé).
#
#     wsl_affected_libraries(<fichier> <variable>)
#
# <fichier> liste les chemins modifiés, relatifs à la racine du dépôt, un par ligne ;
# une ligne « ALL » demande tout. <variable> reçoit « ALL », ou la liste des
# bibliothèques à tester : celles qu'un fichier modifié désigne, puis, par fermeture,
# toutes celles qui en dépendent (DEPENDS de wsl_add_library). Les exemples modifiés
# sont rendus dans <variable>_EXAMPLES.
#
# Correspondance : include/wsl/<nom>.hpp et include/wsl/<nom>/ -> <nom> ;
# tools/<nom>/ -> <nom> ; tests/test_<nom>[_…].cpp ou .hpp -> la plus longue bibliothèque
# préfixe ; tests/compile_fail/<nom>/ -> <nom> ; tests/special_contraction/ -> special ; tests/random_engine_vectors.hpp
# -> random ; tests/package/ : rien (job séparé) ; examples/<nom>.cpp -> cet
# exemple ; chemins de cmake/ci_ignored_paths.txt : rien ; tout autre chemin
# (CMakeLists.txt, cmake/, presets, .github/…) -> ALL.

function(wsl_affected_libraries file variable)
    file(STRINGS ${file} paths)
    file(STRINGS ${PROJECT_SOURCE_DIR}/cmake/ci_ignored_paths.txt ignored)
    get_property(libraries GLOBAL PROPERTY WSL_LIBRARY_NAMES)
    set(changed "")
    set(examples "")
    foreach(path IN LISTS paths)
        string(STRIP "${path}" path)
        if(path STREQUAL "")
            continue()
        endif()
        if(path STREQUAL "ALL")
            set(${variable} ALL PARENT_SCOPE)
            return()
        endif()
        set(skip FALSE)
        foreach(pattern IN LISTS ignored)
            if(path MATCHES "${pattern}")
                set(skip TRUE)
            endif()
        endforeach()
        if(skip OR path MATCHES "^tests/package/")
            continue()
        endif()
        set(library "")
        if(path MATCHES "^include/wsl/([a-z0-9_]+)(\\.hpp|/)")
            set(library ${CMAKE_MATCH_1})
        elseif(path MATCHES "^tools/([a-z0-9_]+)/")
            set(library ${CMAKE_MATCH_1})
        elseif(path MATCHES "^tests/test_([a-z0-9_]+)\\.(cpp|hpp)$")
            set(test ${CMAKE_MATCH_1})
            set(length 0)
            foreach(candidate IN LISTS libraries)
                string(LENGTH ${candidate} candidate_length)
                if((test STREQUAL candidate OR test MATCHES "^${candidate}_")
                   AND candidate_length GREATER length)
                    set(library ${candidate})
                    set(length ${candidate_length})
                endif()
            endforeach()
        elseif(path MATCHES "^tests/compile_fail/([a-z0-9_]+)/")
            set(library ${CMAKE_MATCH_1})
        elseif(path MATCHES "^tests/special_contraction/")
            set(library special)
        elseif(path STREQUAL "tests/random_engine_vectors.hpp")
            set(library random)
        elseif(path MATCHES "^examples/([a-z0-9_]+)\\.cpp$")
            list(APPEND examples ${CMAKE_MATCH_1})
            continue()
        endif()
        if(library STREQUAL "" OR NOT library IN_LIST libraries)
            # Infrastructure ou chemin inconnu : par prudence, tout.
            set(${variable} ALL PARENT_SCOPE)
            return()
        endif()
        list(APPEND changed ${library})
    endforeach()

    # Fermeture : une bibliothèque est touchée si l'une de ses dépendances l'est.
    set(affected ${changed})
    set(grown TRUE)
    while(grown)
        set(grown FALSE)
        foreach(library IN LISTS libraries)
            if(library IN_LIST affected)
                continue()
            endif()
            get_property(dependencies GLOBAL PROPERTY WSL_DEPENDS_${library})
            foreach(dependency IN LISTS dependencies)
                if(dependency IN_LIST affected)
                    list(APPEND affected ${library})
                    set(grown TRUE)
                    break()
                endif()
            endforeach()
        endforeach()
    endwhile()
    list(REMOVE_DUPLICATES affected)
    set(${variable} "${affected}" PARENT_SCOPE)
    set(${variable}_EXAMPLES "${examples}" PARENT_SCOPE)
endfunction()
