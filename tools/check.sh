#!/usr/bin/env bash
# Vérification locale rapide : ne compile et ne teste que les bibliothèques touchées depuis
# main (même règle que la CI d'une PR : cmake/affected.cmake), en release et en FMA par défaut.
#
# Vérifie une **copie figée** de l'état courant (build/check-src, mise à jour par rsync : seuls
# les fichiers modifiés sont copiés, la compilation reste incrémentale) : on peut continuer à
# modifier le dépôt pendant qu'il tourne. Si tout est vert, il enregistre l'arbre vérifié
# (build/check-passed) : le hook d'avant push (tools/git-hooks/pre-push) refuse un push dont le
# code n'est pas cet état.
#
# Usage : tools/check.sh [--all] [--base RÉF] [release|fma|debug]...
#   release (défaut, avec fma) ; fma : Release, -march=native -ffp-contract=fast ;
#   debug : ASan + UBSan ; --base RÉF : comparer à RÉF plutôt qu'à main ; --all : tout tester.
set -euo pipefail
cd "$(dirname "$0")/.."
root=$PWD
base=""
everything=false
configs=()
while (($#)); do
    case "$1" in
        --all) everything=true ;;
        --base) base=$2; shift ;;
        release | fma | debug) configs+=("$1") ;;
        *) echo "usage : tools/check.sh [--all] [--base RÉF] [release|fma|debug]..." >&2; exit 2 ;;
    esac
    shift
done
((${#configs[@]})) || configs=(release fma)
[[ -n "$base" ]] || base=$(git merge-base origin/main HEAD)
mkdir -p build

# L'état vérifié : l'arbre Git de l'arbre de travail (fichiers suivis et nouveaux, hors
# .gitignore), par un index temporaire ; comparable à l'arbre d'un commit.
index=$(mktemp)
trap 'rm -f "$index"' EXIT
GIT_INDEX_FILE=$index git read-tree HEAD
GIT_INDEX_FILE=$index git add -A
tree=$(GIT_INDEX_FILE=$index git write-tree)

# Copie figée : copie les fichiers nouveaux ou modifiés (dates conservées : la compilation
# reste incrémentale), retire ceux qui ont disparu.
snapshot=build/check-src
git ls-files --cached --others --exclude-standard | grep -v '^build/' | python3 -c '
import os, shutil, sys
target = sys.argv[1]
wanted = set(line.rstrip("\n") for line in sys.stdin if os.path.isfile(line.rstrip("\n")))
for path in wanted:
    copy = os.path.join(target, path)
    s = os.stat(path)
    try:
        c = os.stat(copy)
        if c.st_size == s.st_size and int(c.st_mtime) == int(s.st_mtime):
            continue
    except FileNotFoundError:
        os.makedirs(os.path.dirname(copy) or ".", exist_ok=True)
    shutil.copy2(path, copy)
for folder, _, files in os.walk(target):
    for name in files:
        path = os.path.relpath(os.path.join(folder, name), target)
        if path not in wanted:
            os.remove(os.path.join(folder, name))
' "$snapshot"

changed=build/check-changed.txt
if $everything; then
    echo ALL > "$changed"
else
    # Commits de la branche, modifications en cours, fichiers nouveaux.
    { git diff --name-only "$base"; git ls-files --others --exclude-standard; } | sort -u > "$changed"
    if [[ ! -s "$changed" ]]; then
        echo "Rien de modifié depuis main : rien à tester."
        exit 0
    fi
fi
jobs=$(sysctl -n hw.ncpu 2>/dev/null || nproc)
status=0
for config in "${configs[@]}"; do
    dir=build/check-$config
    case "$config" in
        release) preset=release; flags="" ;;
        fma) preset=release; flags="-march=native -ffp-contract=fast" ;;
        debug) preset=debug; flags="" ;;
    esac
    (cd "$snapshot" && cmake --preset "$preset" -S . -B "$root/$dir" \
        -DWSL_AFFECTED_BY="$root/$changed" -DCMAKE_CXX_FLAGS="$flags") > "$dir.configure.log" 2>&1 ||
        { echo "$config : configuration échouée (voir $dir.configure.log)"; status=1; continue; }
    if ! cmake --build "$dir" -j "$jobs" > "$dir.build.log" 2>&1; then
        python3 tools/ci/summary.py job --brief --name "$config" --build-log "$dir.build.log" || true
        echo "  (journal : $dir.build.log)"
        status=1
        continue
    fi
    ctest --test-dir "$dir" -j "$jobs" --no-tests=ignore --output-junit "$root/$dir.junit.xml" > "$dir.test.log" 2>&1 || true
    python3 tools/ci/summary.py job --brief --allow-empty --name "$config" --build-log "$dir.build.log" \
        --junit "$dir.junit.xml" --build-dir "$dir" || status=1
    if grep -q "warning:" "$dir.build.log"; then
        echo "$config : avertissements (voir $dir.build.log)"
        status=1
    fi
done
if ((status == 0)) && [[ " ${configs[*]} " == *" release "* ]]; then
    echo "$tree" > build/check-passed
    echo "État vérifié : $tree"
fi
exit $status
