#pragma once

// wsl::mdarray : tableaux multidimensionnels (conception : docs/mdarray.md).
//
// Lot 1, en cours : positions et déplacements (mdposition, mddisplacement),
// arithmétique entière (div_euclid, rem_euclid), étendues (extents, dextents, dims),
// dispositions (layout_right, layout_left, layout_stride) et vue (mdspan), aux noms
// et à la sémantique de <mdspan> ; conteneur (mdarray, P1684), equal ; axes (axis).

#include <wsl/mdarray/axis.hpp>
#include <wsl/mdarray/container.hpp>
#include <wsl/mdarray/extents.hpp>
#include <wsl/mdarray/integer.hpp>
#include <wsl/mdarray/layout.hpp>
#include <wsl/mdarray/mdspan.hpp>
#include <wsl/mdarray/position.hpp>
