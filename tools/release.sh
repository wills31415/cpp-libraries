#!/usr/bin/env bash
# release.sh — publie la version X.Y.Z, en deux temps (design/interfaces.md, circuit de
# déploiement).
#
#   tools/release.sh prepare X.Y.Z
#       Depuis origin/main : branche release/vX.Y.Z, project(VERSION X.Y.Z), section
#       [Non publié] du CHANGELOG renommée [X.Y.Z] — date (et une nouvelle [Non publié]),
#       liens de comparaison ; commit, push, PR (fusion automatique quand la CI est verte :
#       CMakeLists.txt change, la CI est donc complète).
#
#   tools/release.sh tag X.Y.Z
#       Une fois la PR fusionnée : tag vX.Y.Z annoté sur son commit de fusion, poussé. Le tag
#       ne lance que le bloc Publication de la CI (tools/ci/release_notes.sh), qui crée la
#       release : notes du CHANGELOG et PR incluses, groupées par nature.
set -euo pipefail
info() { printf '\033[36m[release]\033[0m %s\n' "$*"; }
err()  { printf '\033[31m[release] %s\033[0m\n' "$*" >&2; exit 1; }
step="${1:-}"; version="${2:-}"
[[ "$step" =~ ^(prepare|tag)$ && "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] ||
    err "usage : tools/release.sh prepare|tag X.Y.Z"
tag="v$version"; branch="release/$tag"
cd "$(dirname "$0")/.."
git fetch -q origin --tags
git rev-parse -q --verify "refs/tags/$tag" >/dev/null && err "le tag $tag existe déjà"

if [[ "$step" == prepare ]]; then
    [[ -z "$(git status --porcelain)" ]] || err "arbre de travail non propre"
    grep -q '^## \[Non publié\]' CHANGELOG.md || err "pas de section [Non publié] dans CHANGELOG.md"
    git switch -q -c "$branch" origin/main
    sed -i.bak -E "/^project\(/,/\)/s/^( *VERSION )[0-9.]+/\1$version/" CMakeLists.txt && rm CMakeLists.txt.bak
    previous=$(git describe --tags --abbrev=0 origin/main 2>/dev/null || echo "")
    python3 - "$version" "$previous" <<'PY'
import datetime, re, sys
version, previous = sys.argv[1], sys.argv[2]
url = "https://github.com/wills31415/cpp-libraries"
text = open("CHANGELOG.md", encoding="utf-8").read()
today = datetime.date.today().isoformat()
text = text.replace("## [Non publié]", f"## [Non publié]\n\n## [{version}] — {today}", 1)
text = re.sub(r"^\[Non publié\]: .*$", f"[Non publié]: {url}/compare/v{version}...HEAD\n[{version}]: "
              + (f"{url}/compare/{previous}...v{version}" if previous else f"{url}/releases/tag/v{version}"),
              text, count=1, flags=re.M)
open("CHANGELOG.md", "w", encoding="utf-8").write(text)
PY
    git commit -qam "release : version $version"
    git push -q --no-verify -u origin "$branch"   # la PR de version passe la CI complète
    gh pr create --base main --head "$branch" --title "release : version $version" --label interne \
        --body "## Quoi et pourquoi
Version $version : section [$version] du CHANGELOG.

## Changements
- cmake : project(VERSION $version) ;
- CHANGELOG : [Non publié] devient [$version].

## Rupture
Voir la section [$version] du CHANGELOG.

## Vérification
CI complète (CMakeLists.txt modifié). Après la fusion : tools/release.sh tag $version." >/dev/null
    gh pr merge "$branch" --auto --merge >/dev/null
    info "PR de version ouverte, fusion automatique quand la CI est verte : $(gh pr view "$branch" --json url -q .url)"
else
    merge=$(gh pr view "$branch" --json state,mergeCommit -q 'select(.state == "MERGED") | .mergeCommit.oid')
    [[ -n "$merge" ]] || err "la PR $branch n'est pas fusionnée"
    git tag -a "$tag" -m "$tag" "$merge"
    git push -q origin "$tag"
    info "tag $tag poussé sur $(git rev-parse --short "$merge") ; la CI publie la release"
fi
