#!/usr/bin/env bash
# Banc d'essai de wsl::mdarray : un programme par noyau (les mesures se faussent quand
# plusieurs noyaux partagent un programme). Usage : tools/mdarray/benchmark.sh [CXX] [options]
set -euo pipefail
cd "$(dirname "$0")/../.."
cxx=${1:-/usr/bin/clang++}
shift || true
flags=${*:--O3 -march=native}
out=$(mktemp -d)
for kernel in 1 2 3 4 5 6 7 8; do
    "$cxx" -std=c++23 $flags -DNDEBUG -Iinclude -DKERNEL=$kernel tools/mdarray/benchmark.cpp -o "$out/k$kernel"
    "$out/k$kernel"
done
rm -rf "$out"
