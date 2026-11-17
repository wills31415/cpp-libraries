// ENGENDRÉ par tools/random/engine_vectors.py — ne pas modifier à la main.
// Vecteurs de test des moteurs, de sources indépendantes (NumPy, matrices GF(2)).
#pragma once

#include <cstdint>

namespace engine_vectors
{
    // splitmix64, graines 0 et 1234567.
    inline constexpr std::uint64_t splitmix64_seed_0[] = {0xe220a8397b1dcdaf, 0x6e789e6aa1b965f4, 0x06c45d188009454f, 0xf88bb8a8724c81ec, 0x1b39896a51a8749b, 0x53cb9f0c747ea2ea, 0x2c829abe1f4532e1, 0xc584133ac916ab3c};
    inline constexpr std::uint64_t splitmix64_seed_1234567[] = {0x599ed017fb08fc85, 0x2c73f08458540fa5, 0x883ebce5a3f27c77, 0x3fbef740e9177b3f, 0xe3b8346708cb5ecd, 0x6c4f7dbc989944f6, 0x9734aed70f5d5e85, 0x46793dd6f7df31b1};

    // xoshiro256++ : état {1, 2, 3, 4} ; état initialisé par splitmix64(42) ; sauts.
    inline constexpr std::uint64_t xoshiro_state_1234[] = {0x0000000002800001, 0x0000000003800067, 0x000cc00003800067, 0x000cc201994400b2, 0x8012a2019ac433cd, 0x8a69978acdee33ba, 0xc271134733154abd, 0xac2ba09179169e97};
    inline constexpr std::uint64_t xoshiro_seed_42_state[] = {0xbdd732262feb6e95, 0x28efe333b266f103, 0x47526757130f9f52, 0x581ce1ff0e4ae394};
    inline constexpr std::uint64_t xoshiro_seed_42[] = {0xd0764d4f4476689f, 0x519e4174576f3791, 0xfbe07cfb0c24ed8c, 0xb37d9f600cd835b8, 0xcb231c3874846a73, 0x968d9f004e50de7d, 0x201718ff221a3556, 0x9ae94e070ed8cb46};
    inline constexpr std::uint64_t xoshiro_1234_jump_state[] = {0x8c7a153956b5f3d1, 0x701f1a713401d85e, 0x6527f66a65469085, 0x8386b786c4408050};
    inline constexpr std::uint64_t xoshiro_1234_long_jump_state[] = {0x096a8eb71295a400, 0xdbf84991e50f4516, 0x534ee745810d2a0e, 0x31655ca1a2215bf1};

    // pcg64 : état et incrément bruts (ceux de NumPy pour PCG64(0)), puis avance de 10^20.
    inline constexpr std::uint64_t pcg64_state[2] = {0x1aa1b5345996452d, 0x09585eb7a69561e3};   // poids fort, poids faible
    inline constexpr std::uint64_t pcg64_inc[2] = {0x418ddadb3af71a82, 0x588133bc447873a9};   // poids fort, poids faible
    inline constexpr std::uint64_t pcg64_outputs[] = {0xa30febcfd9c2825f, 0x4510bdf882d9d721, 0x0a7d3da94ecde8b8, 0x043b27b61342f01d, 0xd0327a782cde513b, 0xe9aa5979a6401c4e, 0x9b4c7b7180edb27f, 0xbac0495ff8829a45};
    inline constexpr std::uint64_t pcg64_after_advance_1e20[] = {0x3a177933dcb7a373, 0xf44dac2a17c9be6b, 0x5c9d0dee04870149, 0x69221b0fe91085b7};
    // pcg64 initialisé par (seed = 42, stream = 54), à la manière de pcg-cpp.
    inline constexpr std::uint64_t pcg64_seed_42_stream_54[] = {0x86b1da1d72062b68, 0x1304aa46c9853d39, 0xa3670e9e0dd50358, 0xf9090e529a7dae00, 0xc85b9fd837996f2c, 0x606121f8e3919196, 0x7ce1c7ff478354ba, 0xcbc4ac70e541310e};

    // philox4x64 : clé et compteur, sorties dans l'ordre (bloc du compteur, puis compteur + 1).
    inline constexpr std::uint64_t philox_key[] = {0x0123456789abcdef, 0xfedcba9876543210};
    inline constexpr std::uint64_t philox_counter[] = {0x0000000000000005, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000};
    inline constexpr std::uint64_t philox_outputs[] = {0x1419091946cf93c3, 0xd8f38ee79cfe875f, 0x83c953a1723ec495, 0xf3f5cb91950b6a19, 0xd0bba8f1bcf6f692, 0xe3473c643c54e623, 0xeded168e9338e0d9, 0xc20bc8d6143b0f29, 0x1d3f9b580a7a91e2, 0x80b9e338b9d2a202, 0x113bc5c237869222, 0x8494d117fc180028};
    // Compteur dont l'incrément propage une retenue sur deux mots.
    inline constexpr std::uint64_t philox_carry_counter[] = {0xffffffffffffffff, 0xffffffffffffffff, 0x0000000000000000, 0x0000000000000000};
    inline constexpr std::uint64_t philox_carry_outputs[] = {0xe023eda94b604dc0, 0x7645d4a757256b6a, 0x7c617144cf8a0f16, 0x90bea592f697d980, 0x2a6716cc16aa9e72, 0x77a3f007bd213c46, 0x776a367e0881fade, 0xb8f7f9f05d4d1c46};
}
