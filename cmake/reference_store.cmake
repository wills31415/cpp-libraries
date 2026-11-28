# Magasin de références adressé par contenu (mode script) : copie OUTPUT depuis le magasin si
# une famille de même clé y est déjà, sinon la calcule par COMMAND puis l'y range.
#
#     cmake -DNAME=<famille> -DOUTPUT=<fichier> -DSTORE=<dossier> -DSALT=<texte>
#           -DINPUTS=<fichier>|<fichier>… -DCOMMAND=<programme>|<argument>… -P reference_store.cmake
#
# Clé : empreinte SHA-256 du contenu des entrées (générateur, modules) et de SALT (version de
# mpmath), jamais des dates : un checkout ou un rebase ne la change pas.
cmake_minimum_required(VERSION 3.25)

string(REPLACE "|" ";" inputs "${INPUTS}")
string(REPLACE "|" ";" command "${COMMAND}")
set(material "${SALT}")
foreach(input IN LISTS inputs)
    file(SHA256 "${input}" digest)
    string(APPEND material ";${digest}")
endforeach()
string(SHA256 key "${material}")
string(SUBSTRING "${key}" 0 16 key)
set(stored "${STORE}/${NAME}-${key}.hpp")

get_filename_component(directory "${OUTPUT}" DIRECTORY)
file(MAKE_DIRECTORY "${directory}")
if(EXISTS "${stored}")
    file(COPY_FILE "${stored}" "${OUTPUT}")
    message(STATUS "Références ${NAME} : restaurées du magasin (${key})")
    return()
endif()

string(TIMESTAMP start "%s")
execute_process(COMMAND ${command} RESULT_VARIABLE failed)
if(failed OR NOT EXISTS "${OUTPUT}")
    file(REMOVE "${OUTPUT}")
    message(FATAL_ERROR "Références ${NAME} : le générateur a échoué (${failed})")
endif()
string(TIMESTAMP stop "%s")
math(EXPR seconds "${stop} - ${start}")
# Écriture atomique : un autre build peut lire le magasin en même temps.
file(MAKE_DIRECTORY "${STORE}")
string(RANDOM LENGTH 8 suffix)
file(COPY_FILE "${OUTPUT}" "${stored}.${suffix}.tmp")
file(RENAME "${stored}.${suffix}.tmp" "${stored}")
message(STATUS "Références ${NAME} : calculées en ${seconds} s, rangées dans le magasin (${key})")
