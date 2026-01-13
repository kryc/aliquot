#include <iostream>
#include <vector>

#include <gmpxx.h>

#include "factors.hpp"
#include "isprime.hpp"
#include "primefactorcache.hpp"
#include "primefactors.hpp"
#include "primes.hpp"

PrimeFactors
prime_factors(
    const mpz_class& n,
    PrimeFactorCache& cache,
    IsPrime& is_prime
)
{
    // Get the span of prime gaps
    auto gaps = get_prime_gaps();

    PrimeFactors prime_factors;

    std::vector<mpz_class> factors;
    mpz_class prime = 7;
    size_t gap_index = 4;
    mpz_class remainder = n;

    // Divide out 2, 3 and 5 first
    while (remainder % 2 == 0) {
        prime_factors.add_factor(2);
        remainder /= 2;
    }
    while (remainder % 3 == 0) {
        prime_factors.add_factor(3);
        remainder /= 3;
    }
    while (remainder % 5 == 0) {
        prime_factors.add_factor(5);
        remainder /= 5;
    }

    while (remainder > 1 && gap_index < gaps.size()) {
        // Check if the remainder is prime
        if (is_prime.check_small(remainder)) {
            prime_factors.add_factor(remainder);
            return prime_factors;
        }
        // Check if we have the factor in cache
        auto cached_factors = cache.product_exists(remainder.get_ui());
        if (cached_factors.has_value()) {
            auto cached = cached_factors.value();
            prime_factors.update(cached);
            return prime_factors;
        }
        // Check if prime divides remainder
        while (remainder % prime == 0) {
            prime_factors.add_factor(prime);
            // std::cout << "Found factor: " << prime << std::endl;
            remainder /= prime;
        }

        // Get the next VLE-encoded gap
        uint64_t gap = 0;
        uint8_t shift = 0;
        uint8_t byte;
        do {
            byte = gaps[gap_index];
            gap |= static_cast<uint64_t>(byte & 0x7F) << shift;
            shift += 7;
            gap_index++;
        } while ((byte & 0x80) != 0);
        prime += gap;
    }

    if (remainder != 1) {
        // std::cout << "Continuing with prime " << prime << " for remainder " << remainder << std::endl;
    }

    // If remainder is not 1, then we need to continue
    while (remainder > 1) {
        // Check if prime divides remainder
        while (remainder % prime == 0) {
            prime_factors.add_factor(prime);
            remainder /= prime;
        }
        // Get next prime
        mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
    }

    return prime_factors;
}