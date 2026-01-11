#include <iostream>
#include <vector>

#include <gmpxx.h>

#include "factors.hpp"
#include "primefactors.hpp"
#include "primes.hpp"

PrimeFactors
prime_factors(
    const mpz_class& n
)
{
    // Get the span of prime gaps
    auto gaps = get_prime_gaps();
    if (!gaps.has_value()) {
        throw std::runtime_error("Failed to load prime gaps");
    }

    PrimeFactors prime_factors;

    std::vector<mpz_class> factors;
    mpz_class prime = 2;
    size_t gap_index = 1;
    mpz_class remainder = n;

    while (remainder > 1 && gap_index < gaps->size()) {
        // Check if prime divides remainder
        while (remainder % prime == 0) {
            prime_factors.add_factor(prime);
            // std::cout << "Found factor: " << prime << std::endl;
            remainder /= prime;
        }

        // Get the next VLE-encoded gap
        uint64_t gap = 0;
        uint8_t shift = 0;
        while (true) {
            if (gap_index >= gaps->size()) {
                throw std::runtime_error("Ran out of prime gaps");
            }
            uint8_t byte = (*gaps)[gap_index];
            gap |= static_cast<uint64_t>(byte & 0x7F) << shift;
            shift += 7;
            gap_index++;
            if ((byte & 0x80) == 0) {
                break;
            }
        }
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