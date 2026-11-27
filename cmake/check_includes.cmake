# Contrôle des inclusions (design/interfaces.md, décision 31) : chaque #include <wsl/…> d'un
# en-tête doit relever de sa propre bibliothèque, d'un de ses enfants, ou d'une bibliothèque
# dont elle dépend (DEPENDS, transitivement). Sinon le graphe des dépendances est faux, et
# la CI sauterait des tests.
#
#     cmake -DDATA=<build>/libraries.txt -DINCLUDE=<include> -P check_includes.cmake
#
# DATA, écrit par wsl_finalize_libraries() : une ligne par bibliothèque,
#     nom|chemin|dépendances séparées par des virgules|en-têtes (wsl/…) séparés par des virgules

# En mode script, sans version minimale, CMake 3.28 garde l'ancien sens de IN_LIST.
cmake_minimum_required(VERSION 3.28)

file(STRINGS ${DATA} lines)
set(names "")
foreach(line IN LISTS lines)
    # Une expression régulière plutôt que list(GET) : un champ vide (aucune dépendance)
    # disparaîtrait d'une liste sous CMake 3.28.
    if(NOT line MATCHES "^([^|]*)\\|([^|]*)\\|([^|]*)\\|(.*)$")
        message(FATAL_ERROR "${DATA} : ligne illisible « ${line} »")
    endif()
    set(name "${CMAKE_MATCH_1}")
    set(path "${CMAKE_MATCH_2}")
    set(dependencies "${CMAKE_MATCH_3}")
    set(headers "${CMAKE_MATCH_4}")
    string(REPLACE "::" "__" id "${name}")
    string(REPLACE "," ";" dependencies "${dependencies}")
    string(REPLACE "," ";" headers "${headers}")
    list(APPEND names ${name})
    set(path_${id} ${path})
    set(dependencies_${id} ${dependencies})
    set(headers_${id} ${headers})
endforeach()

# Bibliothèque la plus profonde dont relève <chemin> (sans wsl/ ni .hpp).
function(owner_of relative variable)
    set(best "")
    set(best_length -1)
    foreach(name IN LISTS names)
        string(REPLACE "::" "__" id "${name}")
        set(path ${path_${id}})
        string(LENGTH "${path}" length)
        if(relative MATCHES "^${path}(/|$)" AND length GREATER best_length)
            set(best ${name})
            set(best_length ${length})
        endif()
    endforeach()
    set(${variable} "${best}" PARENT_SCOPE)
endfunction()

set(failures 0)
foreach(name IN LISTS names)
    string(REPLACE "::" "__" id "${name}")
    # Permis : la bibliothèque, ses enfants, ses dépendances (et les leurs, et leurs enfants).
    set(allowed ${name})
    set(grown TRUE)
    while(grown)
        set(grown FALSE)
        foreach(member IN LISTS allowed)
            string(REPLACE "::" "__" member_id "${member}")
            foreach(other IN LISTS names)
                if(other IN_LIST allowed)
                    continue()
                endif()
                string(REPLACE "::" "__" other_id "${other}")
                if(other IN_LIST dependencies_${member_id} OR path_${other_id} MATCHES "^${path_${member_id}}/")
                    list(APPEND allowed ${other})
                    set(grown TRUE)
                endif()
            endforeach()
        endforeach()
    endwhile()
    foreach(header IN LISTS headers_${id})
        file(STRINGS ${INCLUDE}/${header} includes REGEX "^[ \t]*#[ \t]*include[ \t]*<wsl/")
        foreach(include IN LISTS includes)
            string(REGEX REPLACE "^[ \t]*#[ \t]*include[ \t]*<wsl/([^>]+)>.*" "\\1" target "${include}")
            string(REGEX REPLACE "\\.hpp$" "" target "${target}")
            owner_of(${target} owner)
            if(owner STREQUAL "")
                message(SEND_ERROR "${header} : <wsl/${target}> ne relève d'aucune bibliothèque")
                math(EXPR failures "${failures} + 1")
            elseif(NOT owner IN_LIST allowed)
                message(SEND_ERROR "${header} (${name}) inclut <wsl/${target}> (${owner}) sans en dépendre : "
                                   "ajouter DEPENDS ${owner} à wsl_add_library(${name})")
                math(EXPR failures "${failures} + 1")
            endif()
        endforeach()
    endforeach()
endforeach()

if(failures GREATER 0)
    message(FATAL_ERROR "${failures} inclusion(s) hors des dépendances déclarées")
endif()
list(LENGTH names count)
message(STATUS "Inclusions conformes aux dépendances (${count} bibliothèques)")
