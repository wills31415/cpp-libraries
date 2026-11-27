# En-têtes sans extension (docs/interfaces.md, décision 30).
#
# Le dépôt ne contient que des .hpp. Pour chaque en-tête public include/wsl/…/X.hpp, on
# engendre un relais sans extension wsl/…/X qui ne fait que l'inclure :
#
#     #include <wsl/random/engines>     →  #include <wsl/random/engines.hpp>
#
# Un fichier et un dossier ne peuvent porter le même nom au même endroit : le parapluie
# wsl/random (fichier) et le dossier wsl/random/ se heurteraient. Un relais dont le nom
# est aussi celui d'un dossier va donc dans une racine d'inclusion à part, une par
# profondeur (umbrella-1 pour wsl/random, umbrella-2 pour wsl/linalg/solvers…) ; les
# autres vont dans la racine « notions ». Le compilateur cherche les racines dans
# l'ordre : umbrella-1, umbrella-2, …, notions. Dans une racine où wsl/random est un
# fichier, <wsl/random/engines> est introuvable et la recherche passe à la suivante
# (vérifié sur GCC 14, Clang 18 et Apple Clang).
#
# Les dossiers detail/ n'ont pas de relais : on ne les inclut qu'entre en-têtes, avec .hpp.
#
#     wsl_generate_relays(<racine des sources> <dossier de sortie> <variable>)
#
# écrit les relais sous <dossier de sortie> (un fichier n'est réécrit que si son contenu
# change, et les relais périmés sont effacés) et met dans <variable> la liste ordonnée des
# racines (chemins relatifs à <dossier de sortie>).

function(wsl_generate_relays source_root output_root roots_variable)
    file(GLOB_RECURSE headers CONFIGURE_DEPENDS RELATIVE ${source_root} ${source_root}/wsl/*.hpp)

    set(expected "")
    set(max_depth 0)
    foreach(header IN LISTS headers)
        if(header MATCHES "(^|/)detail/")
            continue()
        endif()
        string(REGEX REPLACE "\\.hpp$" "" relay ${header})        # wsl/random/engines
        # Profondeur : nombre de composantes après wsl/.
        string(REPLACE "/" ";" parts ${relay})
        list(LENGTH parts length)
        math(EXPR depth "${length} - 1")
        if(IS_DIRECTORY ${source_root}/${relay})
            set(root umbrella-${depth})
            if(depth GREATER max_depth)
                set(max_depth ${depth})
            endif()
        else()
            set(root notions)
        endif()
        set(output ${output_root}/${root}/${relay})
        file(CONFIGURE OUTPUT ${output}
             CONTENT "// -*- C++ -*-\n// Relais engendré par cmake/relays.cmake : ne pas modifier.\n#pragma once\n#include <${header}>\n"
             @ONLY)
        list(APPEND expected ${output})
    endforeach()

    # Relais périmés : un en-tête retiré ou déplacé ne doit pas laisser son relais.
    file(GLOB_RECURSE existing ${output_root}/*)
    foreach(file IN LISTS existing)
        if(NOT file IN_LIST expected)
            file(REMOVE ${file})
        endif()
    endforeach()

    set(roots "")
    if(max_depth GREATER 0)
        foreach(depth RANGE 1 ${max_depth})
            list(APPEND roots umbrella-${depth})
            file(MAKE_DIRECTORY ${output_root}/umbrella-${depth})
        endforeach()
    endif()
    list(APPEND roots notions)
    file(MAKE_DIRECTORY ${output_root}/notions)
    set(${roots_variable} ${roots} PARENT_SCOPE)
endfunction()
