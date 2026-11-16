// ENGENDRÉ par tools/special/generate.py — ne pas modifier à la main.
#pragma once

#include <cstdint>

// Constantes et tables de sin, cos, tan, atan, atan2 (voir trigonometric.hpp).

namespace wsl::special::detail
{
    // π/2 = pio2_1 + … + pio2_5 à 2^-165 près, chaque morceau sur 33 bits.
    inline constexpr double pio2_1 = 0x1.921fb54400000p+0;
    inline constexpr double pio2_2 = 0x1.0b4611a600000p-34;
    inline constexpr double pio2_3 = 0x1.3198a2e000000p-69;
    inline constexpr double pio2_4 = 0x1.b839a25200000p-104;
    inline constexpr double pio2_5 = 0x1.2704453400000p-142;
    inline constexpr double inv_pio2 = 0x1.45f306dc9c883p-1;   // 2/π

    // Constantes en double-double, et arrondis de π/4 et 3π/4.
    inline constexpr double pi_hi = 0x1.921fb54442d18p+1;
    inline constexpr double pi_lo = 0x1.1a62633145c07p-53;
    inline constexpr double pio2_hi = 0x1.921fb54442d18p+0;
    inline constexpr double pio2_lo = 0x1.1a62633145c07p-54;
    inline constexpr double pio4 = 0x1.921fb54442d18p-1;
    inline constexpr double three_pio4 = 0x1.2d97c7f3321d2p+1;

    // Bits de 2/π, 1344 après la virgule, précédés d'un mot nul (Payne–Hanek).
    inline constexpr std::uint64_t two_over_pi_bits[] = {
        0x0000000000000000,
        0xa2f9836e4e441529,
        0xfc2757d1f534ddc0,
        0xdb6295993c439041,
        0xfe5163abdebbc561,
        0xb7246e3a424dd2e0,
        0x06492eea09d1921c,
        0xfe1deb1cb129a73e,
        0xe88235f52ebb4484,
        0xe99c7026b45f7e41,
        0x3991d639835339f4,
        0x9c845f8bbdf9283b,
        0x1ff897ffde05980f,
        0xef2f118b5a0a6d1f,
        0x6d367ecf27cb09b7,
        0x4f463f669e5fea2d,
        0x7527bac7ebe5f17b,
        0x3d0739f78a5292ea,
        0x6bfb5fb11f8d5d08,
        0x56033046fc7b6bab,
        0xf0cfbc209af4361d,
        0xa9e391615ee61b08,
    };

    // Pour c = i / 64, i = 0 … 52 : sin(c) = hi + lo, cos(c) = hi + lo.
    inline constexpr double sincos_table[][4] = {
        {0.0, 0.0, 0x1.0000000000000p+0, 0.0},   // i = 0
        {0x1.fffaaaaeeeed5p-7, -0x1.2ab639a9f0776p-63, 0x1.fff000155549fp-1, 0x1.28a28a03a5ef3p-55},   // i = 1
        {0x1.ffeaaaeeee86fp-6, -0x1.cd406fb224ae2p-60, 0x1.ffc00155527d3p-1, -0x1.3b54492d89b5bp-55},   // i = 2
        {0x1.7fdc01032fba9p-5, -0x1.599bdf46e997ap-59, 0x1.ff7006bfdf99fp-1, -0x1.8b3b560648d5fp-56},   // i = 3
        {0x1.ffaaaeeed4edbp-5, -0x1.2d16d32684b69p-59, 0x1.ff0015549f4d3p-1, 0x1.328387b99426fp-55},   // i = 4
        {0x1.3facb12d1755bp-4, -0x1.921915299468bp-58, 0x1.fe7034129ef6fp-1, -0x1.cbf4337c96f97p-57},   // i = 5
        {0x1.7f701032550e4p-4, 0x1.afc2d1800501ap-60, 0x1.fdc06bf7e6b9bp-1, 0x1.31902b535f8dbp-55},   // i = 6
        {0x1.bf1b78568391dp-4, 0x1.e91841dea4cc8p-58, 0x1.fcf0c800e99b1p-1, 0x1.ea3d786d186acp-57},   // i = 7
        {0x1.feaaeee86ee36p-4, -0x1.afcb2bcc6f03bp-59, 0x1.fc015527d5bd3p-1, 0x1.b68f35094efb8p-55},   // i = 8
        {0x1.1f0d3d7afceafp-3, -0x1.6ef95099769a5p-57, 0x1.faf22263c4bd3p-1, -0x1.52ace133a2769p-58},   // i = 9
        {0x1.3eb312c5d66cbp-3, 0x1.47d666b66cb91p-57, 0x1.f9c340a7cc428p-1, 0x1.c5b6b063b7462p-55},   // i = 10
        {0x1.5e44fcfa126f3p-3, -0x1.6f443063f89b6p-57, 0x1.f874c2e1eecf6p-1, -0x1.c6514e1332b16p-55},   // i = 11
        {0x1.7dc102fbaf2b5p-3, 0x1.5ab50e23c97c3p-59, 0x1.f706bdf9ece1cp-1, -0x1.698c80c36dcb4p-55},   // i = 12
        {0x1.9d252d0cec312p-3, 0x1.9c43d80b1137dp-58, 0x1.f57948cff6797p-1, 0x1.e3a0d3e03b1d4p-57},   // i = 13
        {0x1.bc6f84edc6199p-3, 0x1.9c1a56a7b0cabp-57, 0x1.f3cc7c3b3d16ep-1, -0x1.21a3ad28a3494p-57},   // i = 14
        {0x1.db9e15fb5a5d0p-3, -0x1.32e20d6cc6fc2p-57, 0x1.f20073086649fp-1, 0x1.b940416c1984bp-56},   // i = 15
        {0x1.faaeed4f31577p-3, -0x1.15d88508e32b8p-57, 0x1.f01549f7deea1p-1, 0x1.d3c1e99e5cafdp-55},   // i = 16
        {0x1.0cd00cef36436p-2, -0x1.9fb0a0c93e2b4p-56, 0x1.ee0b1fbc0f11cp-1, -0x1.bfd2380bbc3b1p-59},   // i = 17
        {0x1.1c37d64c6b876p-2, 0x1.46076fe0dcff4p-56, 0x1.ebe214f76efa8p-1, -0x1.02f9f12ba543ep-55},   // i = 18
        {0x1.2b8ddc43eb49fp-2, 0x1.1553899f2d807p-57, 0x1.e99a4c3a7cd83p-1, -0x1.2264b1bc53ce8p-55},   // i = 19
        {0x1.3ad129769d3d8p-2, 0x1.03d550487839ap-63, 0x1.e733ea0193d40p-1, -0x1.6428b3546ce13p-55},   // i = 20
        {0x1.4a00c9b0f3d20p-2, 0x1.823ba6bb08eadp-56, 0x1.e4af14b2a449cp-1, -0x1.68ca02e8a6833p-55},   // i = 21
        {0x1.591bc9fa2f597p-2, 0x1.7c74bac3fe0cbp-57, 0x1.e20bf49acd6c1p-1, -0x1.660aec7ef636bp-58},   // i = 22
        {0x1.682138a38d7f7p-2, -0x1.d889202444aadp-56, 0x1.df4ab3ebd875ep-1, -0x1.e2d8a7e6736c4p-55},   // i = 23
        {0x1.7710255764214p-2, -0x1.6ead7314bb6cep-57, 0x1.dc6b7eb995912p-1, 0x1.4b364776dcd35p-58},   // i = 24
        {0x1.85e7a12826949p-2, 0x1.8a40e9b5face0p-56, 0x1.d96e82f71a9dcp-1, 0x1.ff61bd5d2039dp-55},   // i = 25
        {0x1.94a6be9f546c5p-2, -0x1.69ce13e683f58p-56, 0x1.d653f073e4040p-1, -0x1.76236434bec37p-55},   // i = 26
        {0x1.a34c91cc50ccap-2, -0x1.a310e3b50cecdp-58, 0x1.d31bf8d8d7c06p-1, 0x1.e60dd3089cbddp-56},   // i = 27
        {0x1.b1d8305321617p-2, -0x1.ae242cb99f519p-56, 0x1.cfc6cfa52ad9fp-1, 0x1.8b5b5508f2a0dp-55},   // i = 28
        {0x1.c048b17b140a3p-2, 0x1.19fe6757e9fa7p-57, 0x1.cc54aa2b2972ep-1, 0x1.4ee162ba83a98p-57},   // i = 29
        {0x1.ce9d2e3d4a51fp-2, -0x1.2fc8a12dae298p-57, 0x1.c8c5bf8ce1a84p-1, 0x1.ab3d1a1590123p-56},   // i = 30
        {0x1.dcd4c15329c9ap-2, 0x1.0d4c6e171fd9ap-56, 0x1.c51a48b8b175ep-1, -0x1.1bbb43b9aa880p-57},   // i = 31
        {0x1.eaee8744b05f0p-2, -0x1.789b43c9b027dp-58, 0x1.c1528065b7d50p-1, -0x1.892111312e828p-55},   // i = 32
        {0x1.f8e99e76abc97p-2, 0x1.9d950af2d00a3p-58, 0x1.bd6ea310294f5p-1, 0x1.31bbcc88c109dp-56},   // i = 33
        {0x1.0362939c69955p-1, -0x1.2d8cd78397b01p-55, 0x1.b96eeef58840ep-1, 0x1.45a3cc78fade0p-58},   // i = 34
        {0x1.0a4021e9e1001p-1, -0x1.6f643a13914f6p-55, 0x1.b553a410c104ep-1, 0x1.8ff7947027a15p-58},   // i = 35
        {0x1.110d0c4b69c3bp-1, 0x1.d918998809981p-55, 0x1.b11d04162a4c6p-1, 0x1.1dd561efbc0c2p-56},   // i = 36
        {0x1.17c8e5f2eedb0p-1, 0x1.35e57102e2488p-57, 0x1.accb526f69de5p-1, 0x1.8fb6a8dd6b6ccp-55},   // i = 37
        {0x1.1e7343236574cp-1, 0x1.22a3fa4f41d5ap-56, 0x1.a85ed4373e02dp-1, 0x1.9be06385ec792p-57},   // i = 38
        {0x1.250bb93788bbbp-1, 0x1.ea3d02457bccep-56, 0x1.a3d7d0352bdcfp-1, -0x1.68dbaeca19669p-55},   // i = 39
        {0x1.2b91dea88421ep-1, -0x1.fa371db216ab0p-55, 0x1.9f368ed912f85p-1, -0x1.1d200c5791606p-55},   // i = 40
        {0x1.32054b148bc4fp-1, 0x1.f6b42095a135bp-55, 0x1.9a7b5a36a6514p-1, 0x1.722cfcc9fa7a9p-55},   // i = 41
        {0x1.386597456282bp-1, -0x1.10fada93b07a8p-56, 0x1.95a67e00cb1fdp-1, -0x1.0befda21f862dp-55},   // i = 42
        {0x1.3eb25d36cd53ap-1, -0x1.be570e1570fc0p-58, 0x1.90b84784ddaf7p-1, -0x1.0feb10ab93b87p-56},   // i = 43
        {0x1.44eb381cf386bp-1, -0x1.3ed6c1e6a5505p-55, 0x1.8bb105a5dc900p-1, 0x1.863e03e9474c1p-55},   // i = 44
        {0x1.4b0fc46aab761p-1, 0x1.0da05738cc59cp-61, 0x1.869108d77a6c6p-1, 0x1.338ffe2bfe9ddp-56},   // i = 45
        {0x1.511f9fd7b351cp-1, -0x1.5c0e861c48831p-55, 0x1.8158a31916d5dp-1, -0x1.de8b90b8228dep-57},   // i = 46
        {0x1.571a6966d59b3p-1, 0x1.c843b4d0fb197p-58, 0x1.7c0827f09e54fp-1, -0x1.c73d6d72aee68p-57},   // i = 47
        {0x1.5cffc16bf8f0dp-1, 0x1.96cb370eb578ap-55, 0x1.769fec655211fp-1, -0x1.827d5cf8c68c5p-57},   // i = 48
        {0x1.62cf49921ac79p-1, -0x1.edd9855b6241ap-55, 0x1.712046fa77678p-1, 0x1.425b0a5029c81p-55},   // i = 49
        {0x1.6888a4e134b2fp-1, -0x1.6b7d37644d5e6p-55, 0x1.6b898fa9efb5dp-1, 0x1.15ac786ccf4b2p-56},   // i = 50
        {0x1.6e2b77c40bde1p-1, -0x1.0e729857fad53p-56, 0x1.65dc1fdeb8cbap-1, -0x1.97c1b47337c77p-58},   // i = 51
        {0x1.73b7680dea578p-1, -0x1.2248306dc12a2p-56, 0x1.6018526f563dfp-1, 0x1.46ca5e0e432d0p-55},   // i = 52
    };

    // Pour c = i / 64, i = 0 … 64 : atan(c) = hi + lo.
    inline constexpr double atan_table[][2] = {
        {0.0, 0.0},   // i = 0
        {0x1.fff555bbb729bp-7, -0x1.220c39d4dff50p-61},   // i = 1
        {0x1.ffd55bba97625p-6, -0x1.5ec431444912cp-60},   // i = 2
        {0x1.7fb818430da2ap-5, -0x1.86ef8f794f105p-63},   // i = 3
        {0x1.ff55bb72cfdeap-5, -0x1.c934d86d23f1dp-60},   // i = 4
        {0x1.3f59f0e7c559dp-4, 0x1.ac4ce285df847p-58},   // i = 5
        {0x1.7ee182602f10fp-4, -0x1.cfb654c0c3d98p-58},   // i = 6
        {0x1.be39ebe6f07c3p-4, 0x1.f7b8f29a05987p-58},   // i = 7
        {0x1.fd5ba9aac2f6ep-4, -0x1.cd37686760c17p-59},   // i = 8
        {0x1.1e1fafb043727p-3, -0x1.b485914dacf8cp-59},   // i = 9
        {0x1.3d6eee8c6626cp-3, 0x1.61a3b0ce9281bp-57},   // i = 10
        {0x1.5c9811e3ec26ap-3, -0x1.054ab2c010f3dp-58},   // i = 11
        {0x1.7b97b4bce5b02p-3, 0x1.347b0b4f881cap-58},   // i = 12
        {0x1.9a6a8e96c8626p-3, 0x1.cf601e7b4348ep-59},   // i = 13
        {0x1.b90d7529260a2p-3, 0x1.17b10d2e0e5abp-61},   // i = 14
        {0x1.d77d5df205736p-3, 0x1.c648d1534597ep-57},   // i = 15
        {0x1.f5b75f92c80ddp-3, 0x1.8ab6e3cf7afbdp-57},   // i = 16
        {0x1.09dc597d86362p-2, 0x1.62e47390cb865p-56},   // i = 17
        {0x1.18bf5a30bf178p-2, 0x1.30ca4748b1bf9p-57},   // i = 18
        {0x1.278372057ef46p-2, -0x1.077cdd36dfc81p-56},   // i = 19
        {0x1.362773707ebccp-2, -0x1.963a544b672d8p-57},   // i = 20
        {0x1.44aa436c2af0ap-2, -0x1.5d5e43c55b3bap-56},   // i = 21
        {0x1.530ad9951cd4ap-2, -0x1.2566480884082p-57},   // i = 22
        {0x1.614840309cfe2p-2, -0x1.a725715711f00p-56},   // i = 23
        {0x1.6f61941e4def1p-2, -0x1.c63aae6f6e918p-56},   // i = 24
        {0x1.7d5604b63b3f7p-2, 0x1.69c885c2b249ap-56},   // i = 25
        {0x1.8b24d394a1b25p-2, 0x1.b6d0ba3748fa8p-56},   // i = 26
        {0x1.98cd5454d6b18p-2, 0x1.9e6c988fd0a77p-56},   // i = 27
        {0x1.a64eec3cc23fdp-2, -0x1.24dec1b50b7ffp-56},   // i = 28
        {0x1.b3a911da65c6cp-2, 0x1.ae187b1ca5040p-56},   // i = 29
        {0x1.c0db4c94ec9f0p-2, -0x1.cc1ce70934c34p-56},   // i = 30
        {0x1.cde53432c1351p-2, -0x1.a2cfa4418f1adp-56},   // i = 31
        {0x1.dac670561bb4fp-2, 0x1.a2b7f222f65e2p-56},   // i = 32
        {0x1.e77eb7f175a34p-2, 0x1.0e53dc1bf3435p-56},   // i = 33
        {0x1.f40dd0b541418p-2, -0x1.a3992dc382a23p-57},   // i = 34
        {0x1.0039c73c1a40cp-1, -0x1.b32c949c9d593p-55},   // i = 35
        {0x1.0657e94db30d0p-1, -0x1.d5b495f6349e6p-56},   // i = 36
        {0x1.0c6145b5b43dap-1, 0x1.974fa13b5404fp-58},   // i = 37
        {0x1.1255d9bfbd2a9p-1, -0x1.2bdaee1c0ee35p-58},   // i = 38
        {0x1.1835a88be7c13p-1, 0x1.c621cec00c301p-55},   // i = 39
        {0x1.1e00babdefeb4p-1, -0x1.928df287a668fp-58},   // i = 40
        {0x1.23b71e2cc9e6ap-1, 0x1.c421c9f38224ep-57},   // i = 41
        {0x1.2958e59308e31p-1, -0x1.09e73b0c6c087p-56},   // i = 42
        {0x1.2ee628406cbcap-1, 0x1.c5d5e9ff0cf8dp-55},   // i = 43
        {0x1.345f01cce37bbp-1, 0x1.1021137c71102p-55},   // i = 44
        {0x1.39c391cd4171ap-1, -0x1.2304331d8bf46p-55},   // i = 45
        {0x1.3f13fb89e96f4p-1, 0x1.ecf8b492644f0p-56},   // i = 46
        {0x1.445065b795b56p-1, -0x1.f76d0163f79c8p-56},   // i = 47
        {0x1.4978fa3269ee1p-1, 0x1.2419a87f2a458p-56},   // i = 48
        {0x1.4e8de5bb6ec04p-1, 0x1.4a33dbeb3796cp-55},   // i = 49
        {0x1.538f57b89061fp-1, -0x1.1bb74abda520cp-55},   // i = 50
        {0x1.587d81f732fbbp-1, -0x1.5e5c9d8c5a950p-56},   // i = 51
        {0x1.5d58987169b18p-1, 0x1.0028e4bc5e7cap-57},   // i = 52
        {0x1.6220d115d7b8ep-1, -0x1.2b785350ee8c1p-57},   // i = 53
        {0x1.66d663923e087p-1, -0x1.6ea6febe8bbbap-56},   // i = 54
        {0x1.6b798920b3d99p-1, -0x1.a80386188c50ep-55},   // i = 55
        {0x1.700a7c5784634p-1, -0x1.8c34d25aadef6p-56},   // i = 56
        {0x1.748978fba8e0fp-1, 0x1.7b2a6165884a1p-59},   // i = 57
        {0x1.78f6bbd5d315ep-1, 0x1.406a089803740p-55},   // i = 58
        {0x1.7d528289fa093p-1, 0x1.560821e2f3aa9p-55},   // i = 59
        {0x1.819d0b7158a4dp-1, -0x1.bf76229d3b917p-56},   // i = 60
        {0x1.85d69576cc2c5p-1, 0x1.6b66e7fc8b8c3p-57},   // i = 61
        {0x1.89ff5ff57f1f8p-1, -0x1.55b9a5e177a1bp-55},   // i = 62
        {0x1.8e17aa99cc05ep-1, -0x1.ec182ab042f61p-56},   // i = 63
        {0x1.921fb54442d18p-1, 0x1.1a62633145c07p-55},   // i = 64
    };
}
