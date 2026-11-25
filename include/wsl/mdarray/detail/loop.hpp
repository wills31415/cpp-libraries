#pragma once

// Indications de vectorisation des boucles de wsl::mdarray.

// Itérations indépendantes en mémoire : chaque opérande est disjointe de la destination ou
// lue exactement à la position écrite (vérifié avant, I12). Sans cette indication, le
// compilateur teste le chevauchement de y lu et de y écrit dans y = a x + y, le test échoue
// et la boucle reste scalaire.
#if defined(__clang__)
#    define WSL_MDARRAY_INDEPENDENT_LOOP _Pragma("clang loop vectorize(assume_safety)")
#elif defined(__GNUC__)
#    define WSL_MDARRAY_INDEPENDENT_LOOP _Pragma("GCC ivdep")
#else
#    define WSL_MDARRAY_INDEPENDENT_LOOP
#endif

// Boucle tenue scalaire : Clang vectorise sinon des lectures à pas non unité par des
// lectures dispersées (gather d'AVX2), bien plus lentes ici (indexed.hpp, contiguous_along).
#if defined(__clang__) && defined(__AVX2__)
#    define WSL_MDARRAY_SCALAR_LOOP _Pragma("clang loop vectorize(disable)")
#else
#    define WSL_MDARRAY_SCALAR_LOOP
#endif
