// ENGENDRÉ par tools/special/generate.py — ne pas modifier à la main.
#pragma once

// Constantes et tables de exp et log (voir elementary.hpp pour leur usage).
// Chaque valeur est l'arrondi correct de la valeur exacte ; les paires
// (hi, lo) sont des double-double : hi + lo vaut la valeur à 2^-106 près.

namespace wsl::special::detail
{
    inline constexpr int exp_table_bits = 7;
    inline constexpr double exp_inv_ln2_n = 0x1.71547652b82fep+7;   // N / ln 2
    inline constexpr double exp_ln2_n_hi = 0x1.62e42fefc0000p-8;     // ln 2 / N, 35 bits
    inline constexpr double exp_ln2_n_lo = -0x1.c610ca86c3899p-44;

    // 2^(j/N) pour j = 0 … N-1.
    inline constexpr double exp_table[][2] = {
        {0x1.0000000000000p+0, 0.0},
        {0x1.0163da9fb3335p+0, 0x1.b61299ab8cdb7p-54},
        {0x1.02c9a3e778061p+0, -0x1.19083535b085dp-56},
        {0x1.04315e86e7f85p+0, -0x1.0a31c1977c96ep-54},
        {0x1.059b0d3158574p+0, 0x1.d73e2a475b465p-55},
        {0x1.0706b29ddf6dep+0, -0x1.c91dfe2b13c27p-55},
        {0x1.0874518759bc8p+0, 0x1.186be4bb284ffp-57},
        {0x1.09e3ecac6f383p+0, 0x1.1487818316136p-54},
        {0x1.0b5586cf9890fp+0, 0x1.8a62e4adc610bp-54},
        {0x1.0cc922b7247f7p+0, 0x1.01edc16e24f71p-54},
        {0x1.0e3ec32d3d1a2p+0, 0x1.03a1727c57b53p-59},
        {0x1.0fb66affed31bp+0, -0x1.b9bedc44ebd7bp-57},
        {0x1.11301d0125b51p+0, -0x1.6c51039449b3ap-54},
        {0x1.12abdc06c31ccp+0, -0x1.1b514b36ca5c7p-58},
        {0x1.1429aaea92de0p+0, -0x1.32fbf9af1369ep-54},
        {0x1.15a98c8a58e51p+0, 0x1.2406ab9eeab0ap-55},
        {0x1.172b83c7d517bp+0, -0x1.19041b9d78a76p-55},
        {0x1.18af9388c8deap+0, -0x1.11023d1970f6cp-54},
        {0x1.1a35beb6fcb75p+0, 0x1.e5b4c7b4968e4p-55},
        {0x1.1bbe084045cd4p+0, -0x1.95386352ef607p-54},
        {0x1.1d4873168b9aap+0, 0x1.e016e00a2643cp-54},
        {0x1.1ed5022fcd91dp+0, -0x1.1df98027bb78cp-54},
        {0x1.2063b88628cd6p+0, 0x1.dc775814a8495p-55},
        {0x1.21f49917ddc96p+0, 0x1.2a97e9494a5eep-55},
        {0x1.2387a6e756238p+0, 0x1.9b07eb6c70573p-54},
        {0x1.251ce4fb2a63fp+0, 0x1.ac155bef4f4a4p-55},
        {0x1.26b4565e27cddp+0, 0x1.2bd339940e9d9p-55},
        {0x1.284dfe1f56381p+0, -0x1.a4c3a8c3f0d7ep-54},
        {0x1.29e9df51fdee1p+0, 0x1.612e8afad1255p-55},
        {0x1.2b87fd0dad990p+0, -0x1.10adcd6381aa4p-59},
        {0x1.2d285a6e4030bp+0, 0x1.0024754db41d5p-54},
        {0x1.2ecafa93e2f56p+0, 0x1.1ca0f45d52383p-56},
        {0x1.306fe0a31b715p+0, 0x1.6f46ad23182e4p-55},
        {0x1.32170fc4cd831p+0, 0x1.a9ce78e18047cp-55},
        {0x1.33c08b26416ffp+0, 0x1.32721843659a6p-54},
        {0x1.356c55f929ff1p+0, -0x1.b5cee5c4e4628p-55},
        {0x1.371a7373aa9cbp+0, -0x1.63aeabf42eae2p-54},
        {0x1.38cae6d05d866p+0, -0x1.e958d3c9904bdp-54},
        {0x1.3a7db34e59ff7p+0, -0x1.5e436d661f5e3p-56},
        {0x1.3c32dc313a8e5p+0, -0x1.efff8375d29c3p-54},
        {0x1.3dea64c123422p+0, 0x1.ada0911f09ebcp-55},
        {0x1.3fa4504ac801cp+0, -0x1.7d023f956f9f3p-54},
        {0x1.4160a21f72e2ap+0, -0x1.ef3691c309278p-58},
        {0x1.431f5d950a897p+0, -0x1.1c7dde35f7999p-55},
        {0x1.44e086061892dp+0, 0x1.89b7a04ef80d0p-59},
        {0x1.46a41ed1d0057p+0, 0x1.c944bd1648a76p-54},
        {0x1.486a2b5c13cd0p+0, 0x1.3c1a3b69062f0p-56},
        {0x1.4a32af0d7d3dep+0, 0x1.9cb62f3d1be56p-54},
        {0x1.4bfdad5362a27p+0, 0x1.d4397afec42e2p-56},
        {0x1.4dcb299fddd0dp+0, 0x1.8ecdbbc6a7833p-54},
        {0x1.4f9b2769d2ca7p+0, -0x1.4b309d25957e3p-54},
        {0x1.516daa2cf6642p+0, -0x1.f768569bd93efp-55},
        {0x1.5342b569d4f82p+0, -0x1.07abe1db13cadp-55},
        {0x1.551a4ca5d920fp+0, -0x1.d689cefede59bp-55},
        {0x1.56f4736b527dap+0, 0x1.9bb2c011d93adp-54},
        {0x1.58d12d497c7fdp+0, 0x1.295e15b9a1de8p-55},
        {0x1.5ab07dd485429p+0, 0x1.6324c054647adp-54},
        {0x1.5c9268a5946b7p+0, 0x1.c4b1b816986a2p-60},
        {0x1.5e76f15ad2148p+0, 0x1.ba6f93080e65ep-54},
        {0x1.605e1b976dc09p+0, -0x1.3e2429b56de47p-54},
        {0x1.6247eb03a5585p+0, -0x1.383c17e40b497p-54},
        {0x1.6434634ccc320p+0, -0x1.c483c759d8933p-55},
        {0x1.6623882552225p+0, -0x1.bb60987591c34p-54},
        {0x1.68155d44ca973p+0, 0x1.038ae44f73e65p-57},
        {0x1.6a09e667f3bcdp+0, -0x1.bdd3413b26456p-54},
        {0x1.6c012750bdabfp+0, -0x1.2895667ff0b0dp-56},
        {0x1.6dfb23c651a2fp+0, -0x1.bbe3a683c88abp-57},
        {0x1.6ff7df9519484p+0, -0x1.83c0f25860ef6p-55},
        {0x1.71f75e8ec5f74p+0, -0x1.16e4786887a99p-55},
        {0x1.73f9a48a58174p+0, -0x1.0a8d96c65d53cp-54},
        {0x1.75feb564267c9p+0, -0x1.0245957316dd3p-54},
        {0x1.780694fde5d3fp+0, 0x1.866b80a02162dp-54},
        {0x1.7a11473eb0187p+0, -0x1.41577ee04992fp-55},
        {0x1.7c1ed0130c132p+0, 0x1.f124cd1164dd6p-54},
        {0x1.7e2f336cf4e62p+0, 0x1.05d02ba15797ep-56},
        {0x1.80427543e1a12p+0, -0x1.27c86626d972bp-54},
        {0x1.82589994cce13p+0, -0x1.d4c1dd41532d8p-54},
        {0x1.8471a4623c7adp+0, -0x1.8d684a341cdfbp-55},
        {0x1.868d99b4492edp+0, -0x1.fc6f89bd4f6bap-54},
        {0x1.88ac7d98a6699p+0, 0x1.994c2f37cb53ap-54},
        {0x1.8ace5422aa0dbp+0, 0x1.6e9f156864b27p-54},
        {0x1.8cf3216b5448cp+0, -0x1.0d55e32e9e3aap-56},
        {0x1.8f1ae99157736p+0, 0x1.5cc13a2e3976cp-55},
        {0x1.9145b0b91ffc6p+0, -0x1.dd6792e582524p-54},
        {0x1.93737b0cdc5e5p+0, -0x1.75fc781b57ebcp-57},
        {0x1.95a44cbc8520fp+0, -0x1.64b7c96a5f039p-56},
        {0x1.97d829fde4e50p+0, -0x1.d185b7c1b85d1p-54},
        {0x1.9a0f170ca07bap+0, -0x1.173bd91cee632p-54},
        {0x1.9c49182a3f090p+0, 0x1.c7c46b071f2bep-56},
        {0x1.9e86319e32323p+0, 0x1.824ca78e64c6ep-56},
        {0x1.a0c667b5de565p+0, -0x1.359495d1cd533p-54},
        {0x1.a309bec4a2d33p+0, 0x1.6305c7ddc36abp-54},
        {0x1.a5503b23e255dp+0, -0x1.d2f6edb8d41e1p-54},
        {0x1.a799e1330b358p+0, 0x1.bcb7ecac563c7p-54},
        {0x1.a9e6b5579fdbfp+0, 0x1.0fac90ef7fd31p-54},
        {0x1.ac36bbfd3f37ap+0, -0x1.f9234cae76cd0p-55},
        {0x1.ae89f995ad3adp+0, 0x1.7a1cd345dcc81p-54},
        {0x1.b0e07298db666p+0, -0x1.bdef54c80e425p-54},
        {0x1.b33a2b84f15fbp+0, -0x1.2805e3084d708p-57},
        {0x1.b59728de5593ap+0, -0x1.c71dfbbba6de3p-54},
        {0x1.b7f76f2fb5e47p+0, -0x1.5584f7e54ac3bp-56},
        {0x1.ba5b030a1064ap+0, -0x1.efcd30e54292ep-54},
        {0x1.bcc1e904bc1d2p+0, 0x1.23dd07a2d9e84p-55},
        {0x1.bf2c25bd71e09p+0, -0x1.efdca3f6b9c73p-54},
        {0x1.c199bdd85529cp+0, 0x1.11065895048ddp-55},
        {0x1.c40ab5fffd07ap+0, 0x1.b4537e083c60ap-54},
        {0x1.c67f12e57d14bp+0, 0x1.2884dff483cadp-54},
        {0x1.c8f6d9406e7b5p+0, 0x1.1acbc48805c44p-56},
        {0x1.cb720dcef9069p+0, 0x1.503cbd1e949dbp-56},
        {0x1.cdf0b555dc3fap+0, -0x1.dd83b53829d72p-55},
        {0x1.d072d4a07897cp+0, -0x1.cbc3743797a9cp-54},
        {0x1.d2f87080d89f2p+0, -0x1.d487b719d8578p-54},
        {0x1.d5818dcfba487p+0, 0x1.2ed02d75b3707p-55},
        {0x1.d80e316c98398p+0, -0x1.11ec18beddfe8p-54},
        {0x1.da9e603db3285p+0, 0x1.c2300696db532p-54},
        {0x1.dd321f301b460p+0, 0x1.2da5778f018c3p-54},
        {0x1.dfc97337b9b5fp+0, -0x1.1a5cd4f184b5cp-54},
        {0x1.e264614f5a129p+0, -0x1.7b627817a1496p-54},
        {0x1.e502ee78b3ff6p+0, 0x1.39e8980a9cc8fp-55},
        {0x1.e7a51fbc74c83p+0, 0x1.2d522ca0c8de2p-54},
        {0x1.ea4afa2a490dap+0, -0x1.e9c23179c2893p-54},
        {0x1.ecf482d8e67f1p+0, -0x1.c93f3b411ad8cp-54},
        {0x1.efa1bee615a27p+0, 0x1.dc7f486a4b6b0p-54},
        {0x1.f252b376bba97p+0, 0x1.3a1a5bf0d8e43p-54},
        {0x1.f50765b6e4540p+0, 0x1.9d3e12dd8a18bp-54},
        {0x1.f7bfdad9cbe14p+0, -0x1.dbb12d006350ap-54},
        {0x1.fa7c1819e90d8p+0, 0x1.74853f3a5931ep-55},
        {0x1.fd3c22b8f71f1p+0, 0x1.2eb74966579e7p-57},
    };

    inline constexpr int log_table_bits = 8;
    inline constexpr int log_table_first = 181;
    inline constexpr double log_ln2_hi = 0x1.62e42fefa3800p-1;     // ln 2, 42 bits
    inline constexpr double log_ln2_lo = 0x1.ef35793c76730p-45;

    // Pour c = j / 256, j = 181 … 362 : ln(c) en double-double, et 1/c arrondi.
    inline constexpr double log_table[][3] = {
        {-0x1.630030b3aac49p-2, -0x1.dc18ce51fff99p-57, 0x1.6a13cd1537290p+0},   // j = 181
        {-0x1.5d5bddf595f30p-2, 0x1.6541148cbb8a2p-56, 0x1.6816816816817p+0},   // j = 182
        {-0x1.57bf753c8d1fbp-2, 0x1.0908d15f88b63p-57, 0x1.661ec6a5122f9p+0},   // j = 183
        {-0x1.522ae0738a3d8p-2, 0x1.8f7e9b38a6979p-57, 0x1.642c8590b2164p+0},   // j = 184
        {-0x1.4c9e09e172c3cp-2, 0x1.123615b147a5dp-58, 0x1.623fa77016240p+0},   // j = 185
        {-0x1.4718dc271c41bp-2, -0x1.8fb4c14c56eefp-60, 0x1.6058160581606p+0},   // j = 186
        {-0x1.419b423d5e8c7p-2, -0x1.0dbb243827392p-57, 0x1.5e75bb8d015e7p+0},   // j = 187
        {-0x1.3c25277333184p-2, 0x1.2ad27e50a8ec6p-56, 0x1.5c9882b931057p+0},   // j = 188
        {-0x1.36b6776be1117p-2, 0x1.324f0e883858ep-58, 0x1.5ac056b015ac0p+0},   // j = 189
        {-0x1.314f1e1d35ce4p-2, 0x1.3d69909e5c3dcp-56, 0x1.58ed2308158edp+0},   // j = 190
        {-0x1.2bef07cdc9354p-2, 0x1.82dad7fd86088p-56, 0x1.571ed3c506b3ap+0},   // j = 191
        {-0x1.269621134db92p-2, -0x1.e0efadd9db02bp-56, 0x1.5555555555555p+0},   // j = 192
        {-0x1.214456d0eb8d4p-2, -0x1.f7ae91aeba60ap-57, 0x1.5390948f40febp+0},   // j = 193
        {-0x1.1bf99635a6b95p-2, 0x1.12aeb84249223p-57, 0x1.51d07eae2f815p+0},   // j = 194
        {-0x1.16b5ccbacfb73p-2, -0x1.66fbd28b40935p-56, 0x1.5015015015015p+0},   // j = 195
        {-0x1.1178e8227e47cp-2, 0x1.0e63a5f01c691p-57, 0x1.4e5e0a72f0539p+0},   // j = 196
        {-0x1.0c42d676162e3p-2, -0x1.162c79d5d11eep-58, 0x1.4cab88725af6ep+0},   // j = 197
        {-0x1.07138604d5862p-2, -0x1.cdb16ed4e9138p-56, 0x1.4afd6a052bf5bp+0},   // j = 198
        {-0x1.01eae5626c691p-2, 0x1.18290bd2932e2p-59, 0x1.49539e3b2d067p+0},   // j = 199
        {-0x1.f991c6cb3b379p-3, -0x1.f665066f980a2p-57, 0x1.47ae147ae147bp+0},   // j = 200
        {-0x1.ef5ade4dcffe6p-3, 0x1.08ab2ddc708a0p-58, 0x1.460cbc7f5cf9ap+0},   // j = 201
        {-0x1.e530effe71012p-3, -0x1.2276041f43042p-59, 0x1.446f86562d9fbp+0},   // j = 202
        {-0x1.db13db0d48940p-3, -0x1.aa11d49f96cb9p-58, 0x1.42d6625d51f87p+0},   // j = 203
        {-0x1.d1037f2655e7bp-3, -0x1.60629242471a2p-57, 0x1.4141414141414p+0},   // j = 204
        {-0x1.c6ffbc6f00f71p-3, 0x1.8e58b2c57a4a5p-57, 0x1.3fb013fb013fbp+0},   // j = 205
        {-0x1.bd087383bd8adp-3, -0x1.dd355f6a516d7p-60, 0x1.3e22cbce4a902p+0},   // j = 206
        {-0x1.b31d8575bce3dp-3, 0x1.6353ab386a94dp-57, 0x1.3c995a47babe7p+0},   // j = 207
        {-0x1.a93ed3c8ad9e3p-3, -0x1.bcafa9de97203p-57, 0x1.3b13b13b13b14p+0},   // j = 208
        {-0x1.9f6c407089664p-3, -0x1.35a19605e67efp-59, 0x1.3991c2c187f63p+0},   // j = 209
        {-0x1.95a5adcf7017fp-3, -0x1.142c507fb7a3dp-58, 0x1.3813813813814p+0},   // j = 210
        {-0x1.8beafeb38fe8cp-3, -0x1.55aa8b6997a40p-58, 0x1.3698df3de0748p+0},   // j = 211
        {-0x1.823c16551a3c2p-3, 0x1.1232ce70be781p-57, 0x1.3521cfb2b78c1p+0},   // j = 212
        {-0x1.7898d85444c73p-3, -0x1.ef8f6ebcfb201p-58, 0x1.33ae45b57bcb2p+0},   // j = 213
        {-0x1.6f0128b756abcp-3, 0x1.8de59c21e166cp-57, 0x1.323e34a2b10bfp+0},   // j = 214
        {-0x1.6574ebe8c133ap-3, 0x1.d34f0f4621bedp-60, 0x1.30d190130d190p+0},   // j = 215
        {-0x1.5bf406b543db2p-3, 0x1.1f5b44c0df7e7p-61, 0x1.2f684bda12f68p+0},   // j = 216
        {-0x1.527e5e4a1b58dp-3, 0x1.71a9682395bfdp-61, 0x1.2e025c04b8097p+0},   // j = 217
        {-0x1.4913d8333b561p-3, 0x1.0d5604930f135p-58, 0x1.2c9fb4d812ca0p+0},   // j = 218
        {-0x1.3fb45a59928ccp-3, 0x1.d87e6a354d056p-57, 0x1.2b404ad012b40p+0},   // j = 219
        {-0x1.365fcb0159016p-3, -0x1.7d411a5b944adp-58, 0x1.29e4129e4129ep+0},   // j = 220
        {-0x1.2d1610c86813ap-3, 0x1.499a3f25af95fp-58, 0x1.288b01288b013p+0},   // j = 221
        {-0x1.23d712a49c202p-3, 0x1.6e38161051d69p-57, 0x1.27350b8812735p+0},   // j = 222
        {-0x1.1aa2b7e23f72ap-3, 0x1.c6ef1d9b2ef7ep-59, 0x1.25e22708092f1p+0},   // j = 223
        {-0x1.1178e8227e47cp-3, 0x1.0e63a5f01c691p-58, 0x1.2492492492492p+0},   // j = 224
        {-0x1.08598b59e3a07p-3, 0x1.dd7009902bf32p-57, 0x1.23456789abcdfp+0},   // j = 225
        {-0x1.fe89139dbd566p-4, 0x1.ac9f4215f9393p-58, 0x1.21fb78121fb78p+0},   // j = 226
        {-0x1.ec739830a1120p-4, 0x1.a2bf991780d3fp-59, 0x1.20b470c67c0d9p+0},   // j = 227
        {-0x1.da727638446a2p-4, -0x1.401fa71733019p-58, 0x1.1f7047dc11f70p+0},   // j = 228
        {-0x1.c885801bc4b23p-4, -0x1.a38cb559a6706p-58, 0x1.1e2ef3b3fb874p+0},   // j = 229
        {-0x1.b6ac88dad5b1cp-4, 0x1.0057eed1ca59fp-59, 0x1.1cf06ada2811dp+0},   // j = 230
        {-0x1.a4e7640b1bc38p-4, 0x1.5b5ca203e4259p-58, 0x1.1bb4a4046ed29p+0},   // j = 231
        {-0x1.9335e5d594989p-4, 0x1.478a85704ccb7p-58, 0x1.1a7b9611a7b96p+0},   // j = 232
        {-0x1.8197e2f40e3f0p-4, -0x1.b9f2dffbeed43p-60, 0x1.19453808ca29cp+0},   // j = 233
        {-0x1.700d30aeac0e1p-4, 0x1.72566212cdd05p-61, 0x1.1811811811812p+0},   // j = 234
        {-0x1.5e95a4d9791cbp-4, -0x1.f38745c5c450ap-58, 0x1.16e0689427379p+0},   // j = 235
        {-0x1.4d3115d207eacp-4, -0x1.769f42c7842ccp-58, 0x1.15b1e5f75270dp+0},   // j = 236
        {-0x1.3bdf5a7d1ee64p-4, -0x1.7a976d3b5b45fp-59, 0x1.1485f0e0acd3bp+0},   // j = 237
        {-0x1.2aa04a44717a5p-4, 0x1.d15d38d2fa3f7p-58, 0x1.135c81135c811p+0},   // j = 238
        {-0x1.1973bd1465567p-4, 0x1.7558367a6acf6p-59, 0x1.12358e75d3033p+0},   // j = 239
        {-0x1.08598b59e3a07p-4, 0x1.dd7009902bf32p-58, 0x1.1111111111111p+0},   // j = 240
        {-0x1.eea31c006b87cp-5, 0x1.3e4fc93b7b66cp-59, 0x1.0fef010fef011p+0},   // j = 241
        {-0x1.ccb73cdddb2ccp-5, 0x1.e48fb0500efd4p-59, 0x1.0ecf56be69c90p+0},   // j = 242
        {-0x1.aaef2d0fb10fcp-5, -0x1.a353bb42e0addp-61, 0x1.0db20a88f4696p+0},   // j = 243
        {-0x1.894aa149fb343p-5, -0x1.a8be97660a23dp-60, 0x1.0c9714fbcda3bp+0},   // j = 244
        {-0x1.67c94f2d4bb58p-5, -0x1.0413e6505e603p-59, 0x1.0b7e6ec259dc8p+0},   // j = 245
        {-0x1.466aed42de3eap-5, 0x1.cdd6f7f4a137ep-59, 0x1.0a6810a6810a7p+0},   // j = 246
        {-0x1.252f32f8d183fp-5, 0x1.947f792615916p-59, 0x1.0953f39010954p+0},   // j = 247
        {-0x1.0415d89e74444p-5, -0x1.c05cf1d753622p-59, 0x1.0842108421084p+0},   // j = 248
        {-0x1.c63d2ec14aaf2p-6, 0x1.ce030a686bd86p-60, 0x1.073260a47f7c6p+0},   // j = 249
        {-0x1.8492528c8cabfp-6, 0x1.d192d0619fa67p-60, 0x1.0624dd2f1a9fcp+0},   // j = 250
        {-0x1.432a925980cc1p-6, 0x1.8cdaf39004192p-60, 0x1.05197f7d73404p+0},   // j = 251
        {-0x1.0205658935847p-6, -0x1.27c8e8416e71fp-60, 0x1.0410410410410p+0},   // j = 252
        {-0x1.82448a388a2aap-7, -0x1.04b16137f09a0p-62, 0x1.03091b51f5e1ap+0},   // j = 253
        {-0x1.010157588de71p-7, -0x1.46662d417ced0p-62, 0x1.0204081020408p+0},   // j = 254
        {-0x1.0080559588b35p-8, -0x1.f96638cf63677p-62, 0x1.0101010101010p+0},   // j = 255
        {0.0, 0.0, 0x1.0000000000000p+0},   // j = 256
        {0x1.ff00aa2b10bc0p-9, 0x1.2821ad5a6d353p-63, 0x1.fe01fe01fe020p-1},   // j = 257
        {0x1.fe02a6b106789p-8, -0x1.e44b7e3711ebfp-67, 0x1.fc07f01fc07f0p-1},   // j = 258
        {0x1.7dc475f810a77p-7, -0x1.16d7687d3df21p-62, 0x1.fa11caa01fa12p-1},   // j = 259
        {0x1.fc0a8b0fc03e4p-7, -0x1.83092c59642a1p-62, 0x1.f81f81f81f820p-1},   // j = 260
        {0x1.3cea44346a575p-6, -0x1.0cb5a902b3a1cp-62, 0x1.f6310aca0dbb5p-1},   // j = 261
        {0x1.7b91b07d5b11bp-6, -0x1.5b602ace3a510p-60, 0x1.f44659e4a4271p-1},   // j = 262
        {0x1.b9fc027af9198p-6, -0x1.0ae69229dc868p-64, 0x1.f25f644230ab5p-1},   // j = 263
        {0x1.f829b0e783300p-6, 0x1.33e3f04f1ef23p-60, 0x1.f07c1f07c1f08p-1},   // j = 264
        {0x1.1b0d98923d980p-5, -0x1.e9ae889bac481p-60, 0x1.ee9c7f8458e02p-1},   // j = 265
        {0x1.39e87b9febd60p-5, -0x1.5bfa937f551bbp-59, 0x1.ecc07b301ecc0p-1},   // j = 266
        {0x1.58a5bafc8e4d5p-5, -0x1.ce55c2b4e2b72p-59, 0x1.eae807aba01ebp-1},   // j = 267
        {0x1.77458f632dcfcp-5, 0x1.18d3ca87b9296p-59, 0x1.e9131abf0b767p-1},   // j = 268
        {0x1.95c830ec8e3ebp-5, 0x1.f5a0e80520bf2p-59, 0x1.e741aa59750e4p-1},   // j = 269
        {0x1.b42dd711971bfp-5, -0x1.eb9759c130499p-60, 0x1.e573ac901e574p-1},   // j = 270
        {0x1.d276b8adb0b52p-5, 0x1.1e3c53257fd47p-61, 0x1.e3a9179dc1a73p-1},   // j = 271
        {0x1.f0a30c01162a6p-5, 0x1.85f325c5bbacdp-59, 0x1.e1e1e1e1e1e1ep-1},   // j = 272
        {0x1.075983598e471p-4, 0x1.80da5333c45b8p-59, 0x1.e01e01e01e01ep-1},   // j = 273
        {0x1.16536eea37ae1p-4, -0x1.79da3e8c22cdap-60, 0x1.de5d6e3f8868ap-1},   // j = 274
        {0x1.253f62f0a1417p-4, -0x1.c125963fc4cfdp-62, 0x1.dca01dca01dcap-1},   // j = 275
        {0x1.341d7961bd1d1p-4, -0x1.b599f227becbbp-58, 0x1.dae6076b981dbp-1},   // j = 276
        {0x1.42edcbea646f0p-4, 0x1.ddd4f935996c9p-59, 0x1.d92f2231e7f8ap-1},   // j = 277
        {0x1.51b073f06183fp-4, 0x1.a49e39a1a8be4p-58, 0x1.d77b654b82c34p-1},   // j = 278
        {0x1.60658a93750c4p-4, -0x1.388458ec21b6ap-58, 0x1.d5cac807572b2p-1},   // j = 279
        {0x1.6f0d28ae56b4cp-4, -0x1.906d99184b992p-58, 0x1.d41d41d41d41dp-1},   // j = 280
        {0x1.7da766d7b12cdp-4, -0x1.eeedfcdd94131p-58, 0x1.d272ca3fc5b1ap-1},   // j = 281
        {0x1.8c345d6319b21p-4, -0x1.4a697ab3424a9p-61, 0x1.d0cb58f6ec074p-1},   // j = 282
        {0x1.9ab42462033adp-4, -0x1.2099e1c184e8ep-59, 0x1.cf26e5c44bfc6p-1},   // j = 283
        {0x1.a926d3a4ad563p-4, 0x1.942f48aa70ea9p-58, 0x1.cd85689039b0bp-1},   // j = 284
        {0x1.b78c82bb0eda1p-4, 0x1.0878cf0327e21p-61, 0x1.cbe6d9601cbe7p-1},   // j = 285
        {0x1.c5e548f5bc743p-4, 0x1.5d617ef8161b1p-60, 0x1.ca4b3055ee191p-1},   // j = 286
        {0x1.d4313d66cb35dp-4, 0x1.790dd951d90fap-58, 0x1.c8b265afb8a42p-1},   // j = 287
        {0x1.e27076e2af2e6p-4, -0x1.61578001e0162p-60, 0x1.c71c71c71c71cp-1},   // j = 288
        {0x1.f0a30c01162a6p-4, 0x1.85f325c5bbacdp-58, 0x1.c5894d10d4986p-1},   // j = 289
        {0x1.fec9131dbeabbp-4, -0x1.5746b9981b36cp-58, 0x1.c3f8f01c3f8f0p-1},   // j = 290
        {0x1.0671512ca596ep-3, 0x1.50c647eb86499p-58, 0x1.c26b5392ea01cp-1},   // j = 291
        {0x1.0d77e7cd08e59p-3, 0x1.9a5dc5e9030acp-57, 0x1.c0e070381c0e0p-1},   // j = 292
        {0x1.14785846742acp-3, 0x1.a28813e3a7f07p-57, 0x1.bf583ee868d8bp-1},   // j = 293
        {0x1.1b72ad52f67a0p-3, 0x1.483023472cd74p-58, 0x1.bdd2b899406f7p-1},   // j = 294
        {0x1.2266f190a5acbp-3, 0x1.f547bf1809e88p-57, 0x1.bc4fd65883e7bp-1},   // j = 295
        {0x1.29552f81ff523p-3, 0x1.301771c407dbfp-57, 0x1.bacf914c1bad0p-1},   // j = 296
        {0x1.303d718e47fd3p-3, -0x1.6b9c7d96091fap-63, 0x1.b951e2b18ff23p-1},   // j = 297
        {0x1.371fc201e8f74p-3, 0x1.de6cb62af18a0p-58, 0x1.b7d6c3dda338bp-1},   // j = 298
        {0x1.3dfc2b0ecc62ap-3, -0x1.ab3a8e7d81017p-58, 0x1.b65e2e3beee05p-1},   // j = 299
        {0x1.44d2b6ccb7d1ep-3, 0x1.9f4f6543e1f88p-57, 0x1.b4e81b4e81b4fp-1},   // j = 300
        {0x1.4ba36f39a55e5p-3, 0x1.68981bcc36756p-57, 0x1.b37484ad806cep-1},   // j = 301
        {0x1.526e5e3a1b438p-3, -0x1.746ff8a470d3ap-57, 0x1.b2036406c80d9p-1},   // j = 302
        {0x1.59338d9982086p-3, -0x1.65d22aa8ad7cfp-58, 0x1.b094b31d922a4p-1},   // j = 303
        {0x1.5ff3070a793d4p-3, -0x1.bc60efafc6f6ep-58, 0x1.af286bca1af28p-1},   // j = 304
        {0x1.66acd4272ad51p-3, -0x1.0900e4e1ea8b2p-58, 0x1.adbe87f94905ep-1},   // j = 305
        {0x1.6d60fe719d21dp-3, -0x1.caae268ecd179p-57, 0x1.ac5701ac5701bp-1},   // j = 306
        {0x1.740f8f54037a5p-3, -0x1.b264062a84cdbp-58, 0x1.aaf1d2f87ebfdp-1},   // j = 307
        {0x1.7ab890210d909p-3, 0x1.be36b2d6a0608p-59, 0x1.a98ef606a63bep-1},   // j = 308
        {0x1.815c0a14357ebp-3, -0x1.4be48073a0564p-58, 0x1.a82e65130e159p-1},   // j = 309
        {0x1.87fa06520c911p-3, -0x1.bf7fdbfa08d9ap-57, 0x1.a6d01a6d01a6dp-1},   // j = 310
        {0x1.8e928de886d41p-3, -0x1.569d851a56770p-57, 0x1.a574107688a4ap-1},   // j = 311
        {0x1.9525a9cf456b4p-3, 0x1.d904c1d4e2e26p-57, 0x1.a41a41a41a41ap-1},   // j = 312
        {0x1.9bb362e7dfb83p-3, 0x1.575e31f003e0cp-57, 0x1.a2c2a87c51ca0p-1},   // j = 313
        {0x1.a23bc1fe2b563p-3, 0x1.93711b07a998cp-59, 0x1.a16d3f97a4b02p-1},   // j = 314
        {0x1.a8becfc882f19p-3, -0x1.e8c37918c39ebp-58, 0x1.a01a01a01a01ap-1},   // j = 315
        {0x1.af3c94e80bff3p-3, -0x1.398cff3641985p-58, 0x1.9ec8e951033d9p-1},   // j = 316
        {0x1.b5b519e8fb5a4p-3, 0x1.ba27fdc19e1a0p-57, 0x1.9d79f176b682dp-1},   // j = 317
        {0x1.bc286742d8cd6p-3, 0x1.4fce744870f55p-58, 0x1.9c2d14ee4a102p-1},   // j = 318
        {0x1.c2968558c18c1p-3, -0x1.73dee38a3fb6bp-57, 0x1.9ae24ea5510dap-1},   // j = 319
        {0x1.c8ff7c79a9a22p-3, -0x1.4f689f8434012p-57, 0x1.999999999999ap-1},   // j = 320
        {0x1.cf6354e09c5dcp-3, 0x1.239a07d55b695p-57, 0x1.9852f0d8ec0ffp-1},   // j = 321
        {0x1.d5c216b4fbb91p-3, 0x1.6e443597e4d40p-57, 0x1.970e4f80cb872p-1},   // j = 322
        {0x1.dc1bca0abec7dp-3, 0x1.834c51998b6fcp-57, 0x1.95cbb0be377aep-1},   // j = 323
        {0x1.e27076e2af2e6p-3, -0x1.61578001e0162p-59, 0x1.948b0fcd6e9e0p-1},   // j = 324
        {0x1.e8c0252aa5a60p-3, -0x1.6e03a39bfc89bp-59, 0x1.934c67f9b2ce6p-1},   // j = 325
        {0x1.ef0adcbdc5936p-3, 0x1.48637950dc20dp-57, 0x1.920fb49d0e229p-1},   // j = 326
        {0x1.f550a564b7b37p-3, 0x1.c5f6dfd018c37p-61, 0x1.90d4f120190d5p-1},   // j = 327
        {0x1.fb9186d5e3e2bp-3, -0x1.caaae64f21acbp-57, 0x1.8f9c18f9c18fap-1},   // j = 328
        {0x1.00e6c45ad501dp-2, -0x1.cb9568ff6feadp-57, 0x1.8e6527af1373fp-1},   // j = 329
        {0x1.0402594b4d041p-2, -0x1.28ec217a5022dp-57, 0x1.8d3018d3018d3p-1},   // j = 330
        {0x1.071b85fcd590dp-2, 0x1.d1707f97bde80p-58, 0x1.8bfce8062ff3ap-1},   // j = 331
        {0x1.0a324e27390e3p-2, 0x1.7dcfde8061c03p-56, 0x1.8acb90f6bf3aap-1},   // j = 332
        {0x1.0d46b579ab74bp-2, 0x1.03ec81c3cbd92p-57, 0x1.899c0f601899cp-1},   // j = 333
        {0x1.1058bf9ae4ad5p-2, 0x1.89fa0ab4cb31dp-58, 0x1.886e5f0abb04ap-1},   // j = 334
        {0x1.136870293a8b0p-2, 0x1.7b66298edd24ap-56, 0x1.87427bcc092b9p-1},   // j = 335
        {0x1.1675cababa60ep-2, 0x1.ce63eab883717p-61, 0x1.8618618618618p-1},   // j = 336
        {0x1.1980d2dd4236fp-2, 0x1.9d3d1b0e4d147p-56, 0x1.84f00c2780614p-1},   // j = 337
        {0x1.1c898c16999fbp-2, -0x1.0e5c62aff1c44p-60, 0x1.83c977ab2beddp-1},   // j = 338
        {0x1.1f8ff9e48a2f3p-2, -0x1.c9fdf9a0c4b07p-56, 0x1.82a4a0182a4a0p-1},   // j = 339
        {0x1.22941fbcf7966p-2, -0x1.76f5eb09628afp-56, 0x1.8181818181818p-1},   // j = 340
        {0x1.2596010df763ap-2, -0x1.0f76c57075e9ep-58, 0x1.8060180601806p-1},   // j = 341
        {0x1.2895a13de86a3p-2, 0x1.7ad24c13f040ep-56, 0x1.7f405fd017f40p-1},   // j = 342
        {0x1.2b9303ab89d25p-2, -0x1.896b5fd852ad4p-56, 0x1.7e225515a4f1dp-1},   // j = 343
        {0x1.2e8e2bae11d31p-2, -0x1.8f4cdb95ebdf9p-56, 0x1.7d05f417d05f4p-1},   // j = 344
        {0x1.31871c9544185p-2, -0x1.51acc4c09b379p-60, 0x1.7beb3922e017cp-1},   // j = 345
        {0x1.347dd9a987d55p-2, -0x1.4dd4c580919f8p-57, 0x1.7ad2208e0ecc3p-1},   // j = 346
        {0x1.3772662bfd85bp-2, -0x1.b5629d8117de7p-59, 0x1.79baa6bb6398bp-1},   // j = 347
        {0x1.3a64c556945eap-2, -0x1.c68651945f97cp-57, 0x1.78a4c8178a4c8p-1},   // j = 348
        {0x1.3d54fa5c1f710p-2, -0x1.e3265c6a1c98dp-56, 0x1.77908119ac60dp-1},   // j = 349
        {0x1.404308686a7e4p-2, -0x1.0bcfb6082ce6dp-56, 0x1.767dce434a9b1p-1},   // j = 350
        {0x1.432ef2a04e814p-2, -0x1.29931715ac903p-56, 0x1.756cac201756dp-1},   // j = 351
        {0x1.4618bc21c5ec2p-2, 0x1.f42decdeccf1dp-56, 0x1.745d1745d1746p-1},   // j = 352
        {0x1.49006804009d1p-2, -0x1.9ffc341f177dcp-57, 0x1.734f0c541fe8dp-1},   // j = 353
        {0x1.4be5f957778a1p-2, -0x1.259b35b04813dp-57, 0x1.724287f46debcp-1},   // j = 354
        {0x1.4ec973260026ap-2, -0x1.42a87d977dc5ep-56, 0x1.713786d9c7c09p-1},   // j = 355
        {0x1.51aad872df82dp-2, 0x1.3927ac19f55e3p-59, 0x1.702e05c0b8170p-1},   // j = 356
        {0x1.548a2c3add263p-2, -0x1.819cf7e308ddbp-57, 0x1.6f26016f26017p-1},   // j = 357
        {0x1.5767717455a6cp-2, 0x1.526adb283660cp-56, 0x1.6e1f76b4337c7p-1},   // j = 358
        {0x1.5a42ab0f4cfe2p-2, -0x1.8ebcb7dee9a3dp-56, 0x1.6d1a62681c861p-1},   // j = 359
        {0x1.5d1bdbf5809cap-2, 0x1.4236383dc7fe1p-56, 0x1.6c16c16c16c17p-1},   // j = 360
        {0x1.5ff3070a793d4p-2, -0x1.bc60efafc6f6ep-57, 0x1.6b1490aa31a3dp-1},   // j = 361
        {0x1.62c82f2b9c795p-2, 0x1.7b7af915300e5p-57, 0x1.6a13cd1537290p-1},   // j = 362
    };
}
