#include <iostream>
#include <span>
#include <tuple>

#include "gmpxx.h"

#include "aliquot.hpp"
#include "isprime.hpp"
#include "primefactorcache.hpp"
#include "primefactors.hpp"
#include "primes.hpp"

std::tuple<mpz_class, PrimeFactors>
SumOfDivisors(
    const mpz_class& N,
    PrimeFactorCache<>& Cache,
    const size_t NumThreads
)
{
    // Get prime factors of N
    auto factors = GetPrimeFactors(N, Cache);
    // Cache the factors
    if (Cache.IsOpen()) {
        Cache.Write(factors);
    }
    // Convert the prime factors to a vector of composite factors
    auto composites = factors.GetComposite();
    // Sum the composite factors excluding n itself
    mpz_class sum = 0;
    for (const auto& comp : composites) {
        if (comp != N) {
            sum += comp;
        }
    }
    return {sum, factors};
}

std::tuple<mpz_class, PrimeFactors>
SumOfDivisors(
    const mpz_class& N
)
{
    PrimeFactorCache<> cache;
    return SumOfDivisors(N, cache, std::thread::hardware_concurrency());
}

const bool
DetectLoop(
    const std::span<const mpz_class> Sequence,
    const mpz_class& NextValue
)
{
    for (const auto& value : Sequence.subspan(0, Sequence.size() - 1)) {
        if (value == NextValue) {
            return true;
        }
    }
    return false;
}

std::vector<mpz_class>
AliquotSequence(
    const mpz_class& N,
    const std::string_view CachePath,
    const bool Verbose,
    const size_t NumThreads
)
{
    PrimeFactorCache<> cache(CachePath);
    std::vector<mpz_class> sequence;
    mpz_class current = N;
    size_t index = 0;
    // Output the starting number
    if (Verbose)
        std::cout <<
            std::setw(5) << index++ <<
            " : " <<
            current << std::endl;

    while (true) {
        auto [sum, factors] = SumOfDivisors(current, cache, NumThreads);
        if (sum == 0) {
            break;
        }
        if (Verbose) {
            std::cout <<
                std::setw(5) << index <<
                " : " <<
                sum << " = " <<
                factors.GetString() << std::endl;
        }
        if ((index == 37 && sum != 1471882804) || (index == 52 && sum != 220578719452) || (index == 92 && sum != 8244565422068579772))
        {
            std::cout << "Debug breakpoint reached at index 37 with sum " << sum << std::endl;
            break;
        }
        sequence.push_back(sum);
        if (sum == current || DetectLoop(sequence, sum)) {
            break;
        }
        current = sum;
        index++;
    }
    return sequence;
}