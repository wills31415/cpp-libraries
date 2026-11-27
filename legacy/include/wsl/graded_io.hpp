#pragma once

// En-tête séparé à dessein : <format> et <ostream> sont parmi les plus lourds
// de la bibliothèque standard, et les types eux-mêmes n'ont aucune raison de
// les imposer à qui ne les affiche jamais.

#include <wsl/graded.hpp>

#include <format>
#include <ostream>
#include <string>

// ---------------------------------------------------------------------------
// Forme compacte : les coordonnées de 0 à degree(), entre accolades.
//
//     {1, 2}        pour {1, 2}, {1, 2, 0} et {1, 2, 0, 0}
//     {0, 0, 3}     les zéros de tête sont significatifs et restent affiches
//     {}            l'élément neutre
//
// S'arrêter à degree() n'est pas qu'une économie de caractères : deux valeurs
// équivalentes s'affichent identiquement, ce qui est la propriété qu'on attend
// d'un affichage.
//
// La spécification '#' donne la forme complète, jusqu'à la taille ou au
// plafond, utile en déverminage pour voir la représentation réelle :
//
//     std::format("{}",   jet<int, 4>{1, 2})  ->  {1, 2}
//     std::format("{:#}", jet<int, 4>{1, 2})  ->  {1, 2, 0, 0}
//
// Une seule spécialisation contrainte couvre les deux arités : l'affichage
// ne dépend que de operator[] et de degree(), et deux valeurs égales
// s'affichent donc identiquement quelle que soit leur arité -- ce qui est la
// forme visible de leur identité observationnelle.
// ---------------------------------------------------------------------------
template <wsl::lexicographic X>
struct std::formatter<X, char>
{
    bool full = false;

    constexpr auto parse(std::format_parse_context &context)
    {
        auto position = context.begin();
        if (position != context.end() && *position == '#')
        {
            full = true;
            ++position;
        }
        if (position != context.end() && *position != '}')
            throw std::format_error{"wsl : specification de format invalide"};
        return position;
    }

    template <class Context>
    auto format(X const &operand, Context &context) const
    {
        auto out = context.out();
        *out++ = '{';

        std::size_t const written =
            operand.degree() == wsl::no_support ? 0 : operand.degree() + 1;
        std::size_t const bound = full ? wsl::detail::graded::capacity_of<X> : written;

        for (std::size_t coordinate = 0; coordinate < bound; ++coordinate)
        {
            if (coordinate != 0)
            {
                *out++ = ',';
                *out++ = ' ';
            }
            out = std::format_to(out, "{}", operand[coordinate]);
        }

        *out++ = '}';
        return out;
    }
};

namespace wsl
{
    template <lexicographic X>
    std::string to_string(X const &operand)
    {
        return std::format("{}", operand);
    }

    // Trouve par ADL, donc utilisable sans qualification ni using.
    template <lexicographic X>
    std::ostream &operator<<(std::ostream &stream, X const &operand)
    {
        return stream << std::format("{}", operand);
    }
}

