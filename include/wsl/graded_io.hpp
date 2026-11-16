#ifndef WSL_GRADED_IO_HPP
#define WSL_GRADED_IO_HPP

// En-tete separe a dessein : <format> et <ostream> sont parmi les plus lourds
// de la bibliotheque standard, et les types eux-memes n'ont aucune raison de
// les imposer a qui ne les affiche jamais.

#include <wsl/graded.hpp>

#include <format>
#include <ostream>
#include <string>

// ---------------------------------------------------------------------------
// Forme compacte : les coordonnees de 0 a degree(), entre accolades.
//
//     {1, 2}        pour {1, 2}, {1, 2, 0} et {1, 2, 0, 0}
//     {0, 0, 3}     les zeros de tete sont significatifs et restent affiches
//     {}            l'element neutre
//
// S'arreter a degree() n'est pas qu'une economie de caracteres : deux valeurs
// equivalentes s'affichent identiquement, ce qui est la propriete qu'on attend
// d'un affichage.
//
// La specification '#' donne la forme complete, jusqu'a la taille ou au
// plafond, utile en deverminage pour voir la representation reelle :
//
//     std::format("{}",   jet<int, 4>{1, 2})  ->  {1, 2}
//     std::format("{:#}", jet<int, 4>{1, 2})  ->  {1, 2, 0, 0}
//
// Une seule specialisation contrainte couvre les deux arites : l'affichage
// ne depend que de operator[] et de degree(), et deux valeurs egales
// s'affichent donc identiquement quelle que soit leur arite -- ce qui est la
// forme visible de leur identite observationnelle.
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
        std::size_t const bound = full ? wsl::detail::capacity_of<X> : written;

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

#endif
