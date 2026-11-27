#!/usr/bin/env bash
# Chemins modifiés que la CI doit tester, dans le fichier $1 (un par ligne ; « ALL »
# pour tout). Lu par CMake (WSL_AFFECTED_BY, cmake/affected.cmake).
#
# - push sur une branche : depuis le dernier passage réussi de la CI sur cette branche
#   (un push annulé ou en échec est ainsi couvert par le suivant) ; sinon, tout ;
# - pull request : depuis la base ;
# - tag, lancement manuel, passage planifié : tout.
#
# Variables de GitHub Actions : GITHUB_EVENT_NAME, GITHUB_REF, GITHUB_REF_NAME,
# GITHUB_BASE_REF, GITHUB_REPOSITORY ; GH_TOKEN pour `gh`.
set -euo pipefail
output=$1

all() { echo ALL > "$output"; echo "Tout est testé : $1"; exit 0; }

case "${GITHUB_EVENT_NAME:-}" in
    pull_request)
        base=$(git merge-base "origin/$GITHUB_BASE_REF" HEAD) || all "base introuvable"
        ;;
    push)
        [[ "$GITHUB_REF" == refs/tags/* ]] && all "tag"
        base=$(gh run list --repo "$GITHUB_REPOSITORY" --workflow ci.yml --branch "$GITHUB_REF_NAME" \
                   --event push --status success --limit 1 --json headSha --jq '.[0].headSha' 2>/dev/null || true)
        [[ -n "$base" && "$base" != "null" ]] || all "aucun passage réussi"
        git cat-file -e "$base^{commit}" 2>/dev/null || all "dernier passage réussi absent de l'historique"
        git merge-base --is-ancestor "$base" HEAD || all "dernier passage réussi hors de l'historique"
        ;;
    *)
        all "événement ${GITHUB_EVENT_NAME:-local}"
        ;;
esac

git diff --name-only "$base" HEAD > "$output"
echo "Depuis $base :"
cat "$output"
