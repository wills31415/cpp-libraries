// ENGENDRÉ par tools/elementary_functions/generate.py — ne pas modifier à la main.
#pragma once

#include <wsl/modes.hpp>

// Constantes et tables de digamma (voir gamma.hpp).

namespace wsl::inline WSL_DEFAULTS::special_functions::detail
{
    // x₀, zéro positif de ψ, en trois parties.
    inline constexpr double digamma_root_hi = 0x1.762d86356be3fp+0;
    inline constexpr double digamma_root_lo = 0x1.b86a722197829p-54;
    inline constexpr double digamma_root_lo2 = 0x1.e0d62a6be90c7p-109;

    inline constexpr double euler_gamma_hi = 0x1.2788cfc6fb619p-1;
    inline constexpr double euler_gamma_lo = -0x1.6cb90701fbfabp-58;

    inline constexpr int digamma_root_degree = 24;
    inline constexpr double digamma_root_half_width = 0x1.0000000000000p-3;
    // ψ(x₀ + t) = Σ c_k t^k, |t| <= 1/8 (reste < 2^-66 relativement) : c_1, c_2 en
    // double-double (hi, lo), puis c_3 … c_24.
    inline constexpr double digamma_root_taylor[] = {
        0x1.ef72bc8ee38acp-1,
        -0x1.3879eb97bf58dp-55,
        -0x1.c563b54aa1a35p-2,
        -0x1.c760306906dfep-56,
        0x1.08b4294d50381p-2,
        -0x1.4fc1317257da8p-3,
        0x1.b9a5b6370f3abp-4,
        -0x1.27baba261cc2cp-4,
        0x1.8fce02b239ca7p-5,
        -0x1.0fa7ec36a7d8fp-5,
        0x1.723d6807edcc0p-6,
        -0x1.f970508e1b6a2p-7,
        0x1.5955caaa962f3p-7,
        -0x1.d828079282eb8p-8,
        0x1.42e1acf81d8dcp-8,
        -0x1.b9afc7cee8a14p-9,
        0x1.2e23345f79aafp-9,
        -0x1.9d626f71d1f7ap-10,
        0x1.1acebbd761089p-10,
        -0x1.82f6345c65b35p-11,
        0x1.08bdae1a261d4p-11,
        -0x1.6a3fddea11304p-12,
        0x1.efacab6fb8985p-13,
        -0x1.531f5dc5eb563p-13,
        0x1.d0080f810fab5p-14,
        -0x1.3d7972688af67p-14,
    };

    inline constexpr double digamma_asymptotic_from = 0x1.4000000000000p+3;
    // B_2k / (2k), k = 1 … ; ψ(y) = log y − 1/(2y) − Σ c_k y^(−2k) pour y >= 10
    // (reste < 2^-66 relativement).
    inline constexpr double digamma_asymptotic[] = {
        0x1.5555555555555p-4,
        -0x1.1111111111111p-7,
        0x1.0410410410410p-8,
        -0x1.1111111111111p-8,
        0x1.f07c1f07c1f08p-8,
        -0x1.5995995995996p-6,
        0x1.5555555555555p-4,
        -0x1.c5e5e5e5e5e5ep-2,
        0x1.86e7f9b9fe6e8p+1,
        -0x1.a74ca514ca515p+4,
        0x1.1975cc0ed7304p+8,
        -0x1.c2f0566566566p+11,
        0x1.ac572aaaaaaabp+15,
        -0x1.dc0b1a5cfbe16p+19,
    };
}
