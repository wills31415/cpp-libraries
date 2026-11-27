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
# Correspondance : include/wsl/<chemin> et tests/<chemin> -> la bibliothèque la plus
# profonde dont le chemin relève (wsl_owner, décision 32) ; tools/<nom>/ -> <nom> ;
# tests/package/ : rien (job séparé) ; examples/<nom>.cpp -> cet
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
        if(path MATCHES "^include/wsl/(.+)$")
            wsl_owner(${CMAKE_MATCH_1} library)
        elseif(path MATCHES "^tools/([a-z0-9_]+)/")
            set(library ${CMAKE_MATCH_1})
        elseif(path MATCHES "^tests/(.+)\\.(cpp|hpp)$")
            # Tests rangés comme les espaces de noms (décision 32).
            wsl_owner(${CMAKE_MATCH_1} library)
        endif()
        if(path MATCHES "^examples/([a-z0-9_]+)\\.cpp$")
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
            wsl_library_id(${library} id)
            get_property(dependencies GLOBAL PROPERTY WSL_DEPENDS_${id})
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
