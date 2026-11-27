#!/usr/bin/env bash
# Vérification locale rapide : ne compile et ne teste que les bibliothèques touchées depuis
# le dernier push (même règle que la CI : cmake/affected.cmake), dans des dossiers à part
# (build/check-<configuration>), sans toucher au cache des presets.
#
# Usage : tools/check.sh [--all] [--base RÉF] [release|fma|debug]...
#   configurations : release (par défaut, avec fma), fma (Release, -march=native
#   -ffp-contract=fast), debug (ASan + UBSan) ;
#   --base RÉF : comparer à RÉF plutôt qu'à la branche suivie (origin/main) ;
#   --all : tout tester.
set -euo pipefail
cd "$(dirname "$0")/.."

base="@{upstream}"
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

mkdir -p build
changed=build/check-changed.txt
if $everything; then
    echo ALL > "$changed"
else
    # Commits non poussés, modifications en cours, fichiers nouveaux.
    { git diff --name-only "$base"; git ls-files --others --exclude-standard; } | sort -u > "$changed"
    if [[ ! -s "$changed" ]]; then
        echo "Rien de modifié depuis $base : rien à tester."
        exit 0
    fi
fi

jobs=$(sysctl -n hw.ncpu 2>/dev/null || nproc)
status=0
for config in "${configs[@]}"; do
    dir=build/check-$config
    case "$config" in
        release) cmake --preset release -B "$dir" -DWSL_AFFECTED_BY="$PWD/$changed" -DCMAKE_CXX_FLAGS="" ;;
        fma) cmake --preset release -B "$dir" -DWSL_AFFECTED_BY="$PWD/$changed" \
                 -DCMAKE_CXX_FLAGS="-march=native -ffp-contract=fast" ;;
        debug) cmake --preset debug -B "$dir" -DWSL_AFFECTED_BY="$PWD/$changed" -DCMAKE_CXX_FLAGS="" ;;
    esac > "$dir.configure.log" 2>&1 || { echo "$config : configuration échouée (voir $dir.configure.log)"; status=1; continue; }
    grep -m1 "Bibliothèques testées" "$dir.configure.log" | sed "s/^-- /$config : /" || true
    if ! cmake --build "$dir" -j "$jobs" > "$dir.build.log" 2>&1; then
        echo "$config : compilation échouée (voir $dir.build.log)"
        grep -E "error|warning" "$dir.build.log" | head -20
        status=1
        continue
    fi
    if grep -q "warning:" "$dir.build.log"; then
        echo "$config : avertissements (voir $dir.build.log)"
        status=1
    fi
    if output=$(ctest --test-dir "$dir" -j "$jobs" --output-on-failure --no-tests=ignore 2>&1); then
        echo "$output" | grep -E "tests passed|No tests" | sed "s/^/$config : /" || echo "$config : aucun test"
    else
        echo "$output" | tail -40
        status=1
    fi
done
exit $status
