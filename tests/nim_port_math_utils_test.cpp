#include <cplib/convolution/algorithm_ntt.hpp>
#include <cplib/convolution/bitwise_or_convolution.hpp>
#include <cplib/fps/chirp_z.hpp>
#include <cplib/set_power_series/set_power_series.hpp>
#include <cplib/math/crt.hpp>
#include <cplib/math/online_combination_prefix_sum.hpp>
#include <cplib/utils/count_subarrays_with_sum.hpp>
#include <cplib/utils/digit_dp.hpp>
#include <cplib/utils/implicit_bfs.hpp>
#include <cplib/utils/implicit_dijkstra.hpp>
#include <cplib/utils/parallel_binary_search.hpp>
#include <cplib/utils/game.hpp>
#include <random>
#include <set>

using namespace cplib;
using M = StaticBarrettModint<998244353>;

void testConvolutions(std::mt19937_64 &rng) {
    for (Int n = 0; n <= 8; ++n) {
        Int size = Int(1) << n;
        std::vector<M> a(size), b(size), subset(size), either(size);
        std::vector<Int> plain(size), plainExpected(size);
        for (Int i = 0; i < size; ++i) {
            a[i] = rng() % 100;
            b[i] = rng() % 100;
            plain[i] = rng() % 20;
        }
        for (Int i = 0; i < size; ++i) {
            for (Int j = 0; j < size; ++j) {
                either[i | j] += a[i] * b[j];
                if ((i & j) == 0) {
                    subset[i | j] += a[i] * b[j];
                    plainExpected[i | j] += plain[i] * plain[j];
                }
            }
        }
        assert(subsetConvolution(a, b) == subset);
        assert(subsetConvolution(plain, plain) == plainExpected);
        assert(bitwiseOrConvolution(a, b) == either);
        auto transformed = a;
        subsetZetaTransform(transformed);
        subsetMobiusTransform(transformed);
        assert(transformed == a);
        supersetZetaTransform(transformed);
        supersetMobiusTransform(transformed);
        assert(transformed == a);
        a[0] = 0;
        assert(setLog(setExp(a)) == a);
        a[0] = 3;
        std::vector<M> unit(size);
        unit[0] = 1;
        assert(subsetConvolution(a, setInv(a)) == unit);
        auto square = subsetConvolution(a, a);
        assert(setPow(a, 2) == square);
        assert(setSqrt(square, M(3)) == a);
        auto composed = polynomialCompositeSetPowerSeries(std::vector<M>{2, 3, 5}, a);
        for (Int i = 0; i < size; ++i)
            assert(composed[i] == (i == 0 ? M(2) : M(0)) + a[i] * 3 + square[i] * 5);
        auto projection = setPowerProjection(a, b, n + 4);
        auto power = unit;
        for (M value : projection) {
            M expected = 0;
            for (Int i = 0; i < size; ++i)
                expected += power[i] * b[i];
            assert(value == expected);
            power = subsetConvolution(power, a);
        }
    }
    assert(subsetConvolution(std::vector<M>{}, std::vector<M>{}).empty());
    for (Int n : {0, 1, 2, 15, 70}) {
        std::vector<M> f(n);
        for (auto &x : f)
            x = rng();
        for (M r : {M(0), M(1), M(-1), M(7)}) {
            M a = 11;
            auto values = chirpZ(f, r, 81, a);
            for (Int i = 0; i < Int(values.size()); ++i) {
                M expected = 0;
                for (Int j = n; j-- > 0;)
                    expected = expected * a + f[j];
                assert(values[i] == expected);
                a *= r;
            }
        }
    }
    using Composite = StaticBarrettModint<12>;
    std::vector<Composite> f{1, 2, 3, 4};
    for (Int r = 0; r < 12; ++r) {
        auto values = chirpZ(f, Composite(r), 10, Composite(2));
        Composite x = 2;
        for (auto value : values) {
            Composite expected = 0;
            for (Int j = f.size(); j-- > 0;)
                expected = expected * x + f[j];
            assert(value == expected);
            x *= r;
        }
    }
    AlgorithmNtt first(998244353, 256), second(167772161, 128);
    std::vector<M> a(70), b(80);
    for (auto &x : a)
        x = rng();
    for (auto &x : b)
        x = rng();
    auto left = first.spectrum(a), right = first.spectrum(b);
    auto other = second.spectrum(std::vector<StaticMontgomeryModint<167772161>>{2, 3, 5});
    assert(second.coefficients<StaticMontgomeryModint<167772161>>(other, 0, 3)[2].val() == 5);
    auto product = first.spectrumProduct(left, right);
    auto actual = first.coefficients<M>(product, 0, 149);
    assert(actual == convolution_naive(a, b));
    std::vector<std::uint32_t> accumulated(first.size());
    first.addSpectrumProduct(accumulated, left, right);
    first.addSpectrumProduct(accumulated, left, right);
    auto tail = first.coefficients<M>(accumulated, 128, 21);
    for (Int i = 0; i < 21; ++i)
        assert(tail[i] == actual[i + 128] * 2);
    first.close();
    first.close();
}

void testArithmetic(std::mt19937_64 &rng) {
    for (Int m = 1; m <= 30; ++m) {
        for (Int n = 1; n <= 30; ++n) {
            Int a = Int(rng() % 61) - 30, b = Int(rng() % 61) - 30;
            auto [residue, modulus] = crt({a, b}, {m, n});
            Int expected = -1, period = std::lcm(m, n);
            for (Int x = 0; x < period; ++x)
                if ((x - a) % m == 0 && (x - b) % n == 0) {
                    expected = x;
                    break;
                }
            assert(expected < 0 ? modulus == 0 : residue == expected && modulus == period);
        }
    }
    assert(crt({}, {}) == (std::pair<Int, Int>{0, 1}));
    bool overflow = false;
    try {
        crt({0, 0}, {std::numeric_limits<Int>::max(), 2});
    } catch (const std::overflow_error &) {
        overflow = true;
    }
    assert(overflow);
    for (Int maxN : {0, 1, 2, 10, 101, 500}) {
        auto online = initOnlineCombinationPrefixSum<M>(maxN);
        Combination_Type<M> comb(maxN);
        for (Int n = 0; n <= maxN; ++n) {
            M sum = 0;
            assert(online.prefixSum(n, -1) == M(0));
            for (Int k = 0; k <= n + 1; ++k) {
                sum += comb.ncr(n, k);
                assert(online.prefixSum(n, k) == sum);
            }
            assert(online.rangeSum(n, std::numeric_limits<Int>::min(), std::numeric_limits<Int>::max()) == sum);
        }
    }
    for (Int trial = 0; trial < 300; ++trial) {
        std::vector<Int> a(rng() % 40);
        for (auto &x : a)
            x = Int(rng() % 11) - 5;
        Int target = Int(rng() % 21) - 10, expected = 0;
        for (std::size_t l = 0; l < a.size(); ++l) {
            Int sum = 0;
            for (std::size_t r = l; r < a.size(); ++r) {
                sum += a[r];
                expected += sum == target;
            }
        }
        assert(count_subarrays_with_sum(a, target) == expected);
    }
    assert(count_subarrays_with_sum(std::vector<UInt>{0, 0, 2}, UInt(0)) == 3);
    overflow = false;
    try {
        count_subarrays_with_sum(std::vector<Int>{std::numeric_limits<Int>::max(), 1}, Int(0));
    } catch (const std::overflow_error &) {
        overflow = true;
    }
    assert(overflow);
}

void testUtilities() {
    for (Int upper : {0, 1, 9, 31, 99, 314, 1200}) {
        auto counts = digitDp(upper, Int(0), [](Int, Int sum, Int d) { return (sum + d) % 7; }, Int(0), Int(1));
        std::unordered_map<Int, Int> expected;
        for (Int x = 0; x <= upper; ++x) {
            Int sum = 0;
            for (Int t = x; t; t /= 10)
                sum += t % 10;
            ++expected[sum % 7];
        }
        assert(counts == expected);
        auto initial = [](Int first, Int digits) -> std::optional<Int> {
            if (first == 4)
                return std::nullopt;
            return digits;
        };
        auto transition = [](Int, Int length, Int d) -> std::optional<Int> {
            return d == 4 ? std::nullopt : std::optional<Int>(length);
        };
        auto filtered = digitDp(std::to_string(upper), initial, transition, Int(0), Int(1));
        expected.clear();
        for (Int x = 0; x <= upper; ++x) {
            auto s = std::to_string(x);
            if (s.find('4') == std::string::npos)
                ++expected[s.size()];
        }
        assert(filtered == expected);
    }
    auto adjacent = [](Int x) {
        std::vector<Int> out;
        if (x > 0)
            out.push_back(x - 1);
        if (x < 20)
            out.push_back(x + 1);
        return out;
    };
    std::vector<Int> starts{2, 2, 15};
    auto result = restore_implicit_bfs(starts, adjacent);
    for (Int x = 0; x <= 20; ++x)
        assert(result.costs.at(x) == std::min(std::abs(x - 2), std::abs(x - 15)));
    auto path = shortest_path_implicit_bfs(starts, std::vector<Int>{5, 19}, adjacent);
    assert(path.cost == 3 && path.path.front() == 2 && path.path.back() == 5);
    assert(implicit_bfs_until(Int(2), adjacent, [](Int x) { return x == 2; }) == 0);
    auto weighted = [&](Int x) {
        std::vector<std::pair<Int, Int>> out;
        for (Int y : adjacent(x))
            out.emplace_back(y, y > x ? 2 : 1);
        return out;
    };
    auto distances = implicit_dijkstra(starts, weighted);
    for (Int x = 0; x <= 20; ++x) {
        auto distance = [x](Int s) { return x >= s ? 2 * (x - s) : s - x; };
        assert(distances.at(x) == std::min(distance(2), distance(15)));
    }
    assert(implicit_bfs_until(starts,
        [](auto x) { return x < 20 ? std::vector<Int>{x + 1} : std::vector<Int>{}; },
        [](auto x) { return x == 5; }) == 3);
    assert(implicit_dijkstra_until(starts,
        [](auto x) { return x < 20 ? std::vector<std::pair<Int, Int>>{{x + 1, 2}} : std::vector<std::pair<Int, Int>>{}; },
        [](auto x) { return x == 5; }) == 6);
    auto shortest = shortest_path_implicit_dijkstra(starts, Int(10), weighted);
    assert(shortest.cost == 5 && shortest.path.front() == 15);
    for (Int m : {0, 1, 30}) {
        Int state = 0;
        auto found = parallelBinarySearch(m, m + 3, [&]() { state = 0; },
            [&](Int index) { assert(index == state); ++state; }, [&](Int q) { return state >= q; });
        for (Int q = 0; q < m + 3; ++q)
            assert(found[q] == std::min(q, m + 1));
    }
    auto next = [](Int x) { return x == 0 ? std::vector<Int>{} : std::vector<Int>{x - 1}; };
    auto terminal = [](Int x) { assert(x == 0); return true; };
    assert(can_win(Int(0), next, false, terminal));
    assert(!can_win(Int(1), next, false, terminal));
    assert(optimal_play(Int(2), next, false, terminal).is_win);
}
int main() {
    std::mt19937_64 rng(20261113);
    testConvolutions(rng);
    testArithmetic(rng);
    testUtilities();
}
