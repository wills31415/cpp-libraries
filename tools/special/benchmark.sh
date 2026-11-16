#!/usr/bin/env bash
# benchmark.sh — coût de wsl::special contre <cmath>, en nanosecondes par appel.
#
# Usage : tools/special/benchmark.sh        (compilateur : $CXX, sinon /usr/bin/clang++)
#
# Un programme par fonction mesurée (voir benchmark.cpp), minimum de 9 mesures.
# Chiffres informatifs : ils dépendent de la machine et de sa charge.

set -euo pipefail
cd "$(dirname "$0")/../.."
CXX="${CXX:-/usr/bin/clang++}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

measure() {   # fonction type bas haut
    "$CXX" -std=c++23 -O3 -DNDEBUG -Iinclude "-DFUNCTION=$1" "-DTYPE=$2" "-DLOW=$3" "-DHIGH=$4" \
        tools/special/benchmark.cpp -o "$work/b"
    "$work/b" | cut -d' ' -f1
}

printf '%-12s %10s %10s %8s\n' "" "wsl" "<cmath>" "rapport"
while read -r name fn type low high; do
    ours="$(measure "wsl::special::$fn" "$type" "$low" "$high")"
    theirs="$(measure "std::$fn" "$type" "$low" "$high")"
    printf '%-12s %7s ns %7s ns %7s\n' "$name" "$ours" "$theirs" "×$(echo "scale=2; $ours / $theirs" | bc)"
done <<'TABLE'
exp_double exp double -700 700
log_double log double 1e-6 1e6
exp_float exp float -80 80
log_float log float 1e-6 1e6
TABLE
