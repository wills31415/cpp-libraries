#!/usr/bin/env bash
# release.sh — publie la version X.Y.Z de cpp-libraries.
#
# Usage : tools/release.sh X.Y.Z
#
# Préalable (à faire dans un commit poussé sur main) :
#   - project(... VERSION X.Y.Z) dans CMakeLists.txt ;
#   - une section « ## [X.Y.Z] — AAAA-MM-JJ » dans CHANGELOG.md.
#
# Le script vérifie tout cela, puis que les tests passent en local et que la CI
# est verte sur ce commit, et seulement alors crée le tag vX.Y.Z et la release
# GitHub (notes tirées du CHANGELOG). Les projets le voient aussitôt : cpp-init
# fige les nouveaux projets sur cette version.

set -euo pipefail

info() { printf '\033[36m[release]\033[0m %s\n' "$*"; }
err()  { printf '\033[31m[release] %s\033[0m\n' "$*" >&2; exit 1; }

version="${1:-}"
[[ "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] || err "usage : tools/release.sh X.Y.Z"
tag="v$version"

cd "$(dirname "$0")/.."

[[ "$(git branch --show-current)" == main ]] || err "à lancer depuis main"
[[ -z "$(git status --porcelain)" ]] || err "arbre de travail non propre"
git fetch -q origin
[[ "$(git rev-parse HEAD)" == "$(git rev-parse origin/main)" ]] || err "main n'est pas à jour avec origin/main (pousser d'abord)"
git rev-parse -q --verify "refs/tags/$tag" >/dev/null && err "le tag $tag existe déjà"

cmake_version="$(sed -n 's/^ *VERSION \([0-9.]*\)$/\1/p' CMakeLists.txt | head -1)"
[[ "$cmake_version" == "$version" ]] || err "CMakeLists.txt annonce la version $cmake_version, pas $version"

notes="$(awk -v v="$version" '
    $0 ~ "^## \\[" v "\\]" {inside = 1; next}
    inside && /^## \[/ {exit}
    inside && /^\[[^]]*\]: / {exit}
    inside {print}' CHANGELOG.md)"
[[ -n "${notes//[[:space:]]/}" ]] || err "pas de section [$version] dans CHANGELOG.md"

info "tests locaux (Release)…"
cmake --preset release >/dev/null
cmake --build --preset release -j >/dev/null
ctest --preset release >/dev/null || err "les tests échouent en local"

info "CI sur $(git rev-parse --short HEAD)…"
run="$(gh run list --workflow CI --commit "$(git rev-parse HEAD)" --limit 1 --json databaseId -q '.[0].databaseId')"
[[ -n "$run" ]] || err "aucune exécution de la CI pour ce commit"
gh run watch "$run" --exit-status >/dev/null || err "la CI n'est pas verte (gh run view $run)"

git tag -a "$tag" -m "$tag"
git push -q origin "$tag"
gh release create "$tag" --title "$tag" --notes "$notes" >/dev/null
info "publié : $tag ($(gh release view "$tag" --json url -q .url))"
