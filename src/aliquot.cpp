#include <iostream>
#include <span>

#include "gmpxx.h"

#include "aliquot.hpp"
#include "isprime.hpp"
#include "primefactorcache.hpp"
#include "primefactors.hpp"
#include "primes.hpp"

mpz_class
sum_of_divisors(
    const mpz_class& n,
    PrimeFactorCache& cache,
    const size_t num_threads
)
{
    // Get prime factors of n
    auto factors = prime_factors(n, cache);
    // Convert the prime factors to a vector of composite factors
    auto composites = factors.get_composite();
    // Sum the composite factors excluding n itself
    mpz_class sum = 0;
    for (const auto& comp : composites) {
        if (comp != n) {
            sum += comp;
        }
    }
    return sum;
}

mpz_class
sum_of_divisors(
    const mpz_class& n
)
{
    PrimeFactorCache cache;
    return sum_of_divisors(n, cache, std::thread::hardware_concurrency());
}

const bool
detect_loop(
    const std::span<const mpz_class> sequence,
    const mpz_class& next_value
)
{
    for (const auto& value : sequence.subspan(0, sequence.size() - 1)) {
        if (value == next_value) {
            return true;
        }
    }
    return false;
}

std::vector<mpz_class>
aliquot_sequence(
    const mpz_class& n,
    const std::string_view cache_path,
    const bool verbose,
    const size_t num_threads
)
{
    PrimeFactorCache cache(cache_path);
    std::vector<mpz_class> sequence;
    mpz_class current = n;
    size_t index = 0;
    while (true) {
        mpz_class sum = sum_of_divisors(current, cache, num_threads);
        if (sum == 0) {
            break;
        }
        if (verbose)
            std::cout << index << ": " << sum << std::endl;
        sequence.push_back(sum);
        if (sum == current || detect_loop(sequence, sum)) {
            break;
        }
        current = sum;
        index++;
    }
    return sequence;
}