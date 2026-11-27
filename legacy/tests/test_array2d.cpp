#include <wsl/array2d.hpp>
#include <cassert>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using wsl::array2d, wsl::row_major, wsl::column_major, wsl::strided_view;

using RowMat = array2d<int, 2, 3>;
using ColMat = array2d<int, 2, 3, column_major>;

// ---- Propriétés de type -----------------------------------------------------
static_assert(std::same_as<RowMat, array2d<int, 2, 3, row_major>>);
static_assert(std::is_aggregate_v<RowMat>);
static_assert(std::is_trivially_copyable_v<RowMat>);
static_assert(sizeof(RowMat) == 6 * sizeof(int));
static_assert(std::ranges::contiguous_range<RowMat>);
static_assert(std::ranges::sized_range<RowMat>);

static_assert(std::ranges::random_access_range<strided_view<int>>);
static_assert(std::ranges::sized_range<strided_view<int>>);
static_assert(std::ranges::view<strided_view<int>>);
static_assert(std::ranges::borrowed_range<strided_view<int>>);
static_assert(std::random_access_iterator<strided_view<const int>::iterator>);
static_assert(std::sortable<strided_view<int>::iterator>);

// Ligne majeure : lignes contiguës (span), colonnes à pas (strided_view). Et l'inverse.
static_assert(std::same_as<decltype(std::declval<RowMat&>().row(0)), std::span<int, 3>>);
static_assert(std::same_as<decltype(std::declval<const RowMat&>().row(0)), std::span<const int, 3>>);
static_assert(std::same_as<decltype(std::declval<RowMat&>().column(0)), strided_view<int>>);
static_assert(std::same_as<decltype(std::declval<ColMat&>().column(0)), std::span<int, 2>>);
static_assert(std::same_as<decltype(std::declval<ColMat&>().row(0)), strided_view<int>>);
static_assert(std::ranges::random_access_range<decltype(std::declval<RowMat&>().rows())>);

// Les vues refusent les temporaires
template <typename A>
concept row_on_temporary = requires { std::declval<A&&>().row(0); };
static_assert(!row_on_temporary<RowMat>);

// ---- Tests à la compilation -------------------------------------------------
template <typename Mat>
consteval bool scenario() {
    auto a = Mat::from_rows({{1, 2, 3}, {4, 5, 6}});
    if (a[1, 2] != 6 || a[0, 1] != 2) return false;
    // somme des colonnes
    int column_sums[3]{};
    int j = 0;
    for (auto col : a.columns()) { for (int x : col) column_sums[j] += x; ++j; }
    if (column_sums[0] != 5 || column_sums[1] != 7 || column_sums[2] != 9) return false;
    // somme des lignes
    int row_sums[2]{};
    int i = 0;
    for (auto row : a.rows()) { for (int x : row) row_sums[i] += x; ++i; }
    if (row_sums[0] != 6 || row_sums[1] != 15) return false;
    // parcours logique (ligne par ligne) quelle que soit la disposition
    int expected = 1;
    for (int x : a.rows() | std::views::join) if (x != expected++) return false;
    // transposée
    auto t = a.transposed();
    if (t[2, 1] != 6 || t.extent(0) != 3 || t.extent(1) != 2 || t.rank() != 2) return false;
    if (t.transposed() != a) return false;
    // tri d'une colonne à travers la vue
    a[0, 1] = 9;
    std::ranges::sort(a.column(1));
    if (a[0, 1] != 5 || a[1, 1] != 9) return false;
    // modification à travers une ligne
    for (int& x : a.row(0)) x = 0;
    if (a[0, 0] != 0 || a[0, 2] != 0 || a[1, 0] != 4) return false;
    return true;
}
static_assert(scenario<RowMat>());
static_assert(scenario<ColMat>());

// Ordre mémoire différent selon la disposition
static_assert(RowMat::from_rows({{1, 2, 3}, {4, 5, 6}}).elems_[1] == 2);
static_assert(ColMat::from_rows({{1, 2, 3}, {4, 5, 6}}).elems_[1] == 4);

// Initialisation d'agrégat directe (ordre mémoire)
constexpr RowMat literal{1, 2, 3, 4, 5, 6};
static_assert(literal[1, 0] == 4);

int main() {
    auto a = RowMat::from_rows({{1, 2, 3}, {4, 5, 6}});

    // at vérifie chaque indice ; parenthèses doubles : la virgule de a[i, j]
    // serait sinon vue comme un séparateur d'arguments de la macro assert
    bool thrown = false;
    try { (void)a.at(0, 5); } catch (const std::out_of_range&) { thrown = true; }
    assert(thrown);
    assert((a.at(1, 1) == 5));

    // algorithmes à plat
    assert(std::accumulate(a.begin(), a.end(), 0) == 21);
    assert(std::ranges::max(a) == 6);

    // strided_view : view_interface fournit operator[], front, back
    auto last_column = a.column(2);
    assert(last_column[1] == 6 && last_column.front() == 3 && last_column.back() == 6);
    assert(last_column.size() == 2);

    // types non triviaux
    array2d<std::string, 2, 2> labels{};
    labels[0, 1] = "haut-droite";
    auto copy = labels;
    assert((copy == labels && copy[0, 1] == "haut-droite"));
    labels.fill("x");
    swap(labels, copy);
    assert((labels[0, 1] == "haut-droite" && copy[1, 1] == "x"));

    // matérialiser une colonne
    auto middle = a.column(1) | std::ranges::to<std::vector>();
    assert((middle == std::vector{2, 5}));

#ifdef __cpp_lib_mdspan
    auto view = a.to_mdspan();
    static_assert(decltype(view)::rank() == 2 && decltype(view)::static_extent(1) == 3);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 3; ++j) assert((view[i, j] == a[i, j]));
    auto col_major = ColMat::from_rows({{1, 2, 3}, {4, 5, 6}});
    const auto& col_major_ref = col_major;
    auto const_view = col_major_ref.to_mdspan();
    static_assert(std::is_const_v<typename decltype(const_view)::element_type>);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 3; ++j) assert((const_view[i, j] == col_major[i, j]));
    view[0, 0] = 42;
    assert((a[0, 0] == 42));
    std::cout << "mdspan testé\n";
#else
    std::cout << "mdspan indisponible, conversion non testée\n";
#endif
    std::cout << "OK\n";
}