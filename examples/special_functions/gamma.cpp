// Fonctions spéciales : gamma, fonction d'erreur, bêta incomplète, chacune dans les deux modes,
// comme les fonctions élémentaires.

#include <wsl/special_functions>

#include <cstdio>

namespace sf = wsl::special_functions;

int main()
{
    std::printf("Γ(4.5)          = %.17g\n", sf::tgamma(4.5));
    std::printf("log Γ(100)      = %.17g\n", sf::lgamma(100.0));
    std::printf("erfc(3)         = %.17g\n", sf::erfc(3.0, wsl::portable));
    std::printf("I_0.3(2, 5)     = %.17g\n", sf::ibeta(2.0, 5.0, 0.3));
    double const p = sf::gamma_p(2.5, 1.0);
    std::printf("P(2.5, 1)       = %.17g, réciproque %.17g\n", p, sf::gamma_p_inv(2.5, p));
    return sf::tgamma(5.0) == 24.0 ? 0 : 1;
}
