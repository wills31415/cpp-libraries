#!/usr/bin/env bash
# Publication d'une version (bloc Publication de la CI) : vérifie le tag vX.Y.Z, puis écrit
# sur la sortie la section [X.Y.Z] du CHANGELOG, qui ouvre les notes de version (GitHub y
# ajoute la liste des PR incluses, groupées par nature : .github/release.yml).
#
# Refuse : un tag mal formé ; une version différente de project(VERSION) ; une section
# [X.Y.Z] absente du CHANGELOG ; un commit hors de main. La PR de version a changé
# CMakeLists.txt : elle a donc passé la CI complète.
set -euo pipefail
tag=$1
[[ "$tag" =~ ^v([0-9]+\.[0-9]+\.[0-9]+)$ ]] || { echo "tag mal formé : $tag (attendu vX.Y.Z)" >&2; exit 1; }
version=${BASH_REMATCH[1]}
declared=$(sed -n '/^project(/,/)/s/^ *VERSION \([0-9.]*\).*/\1/p' CMakeLists.txt | head -1)
[[ "$declared" == "$version" ]] || { echo "project(VERSION $declared) ne correspond pas au tag $tag" >&2; exit 1; }
git fetch --quiet origin main
git merge-base --is-ancestor HEAD origin/main || { echo "le commit de $tag n'est pas sur main" >&2; exit 1; }
section=$(awk -v v="$version" '
    $0 ~ "^## \\[" v "\\]" { inside = 1; next }
    inside && /^## \[/ { exit }
    inside { print }' CHANGELOG.md)
[[ -n "${section//[[:space:]]/}" ]] || { echo "section [$version] absente ou vide dans CHANGELOG.md" >&2; exit 1; }
printf '%s\n' "$section"
