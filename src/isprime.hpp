#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include <gmpxx.h>

#include "primes.hpp"

// A class to check whether a number is prime using a bitset for small primes
// And falling back to GMP for larger numbers
class IsPrime {
public:
    IsPrime(void) {
        Build(GetPrimeGaps());
    }
    IsPrime(std::span<const uint8_t> PrimeGaps) {
        Build(PrimeGaps);
    }

    template <typename T>
    const bool
    CheckSmall(
        const T& N
    ) const {
        if (N < 2) {
            return false;
        }
        if (N == 2) {
            return true;
        }
        if constexpr (std::is_same_v<T, mpz_class>) {
            if ((N.get_ui() & 1) == 0) {
                return false;  // Even numbers > 2 are not prime
            }
        } else {
            if ((N & 1) == 0) {
                return false;  // Even numbers > 2 are not prime
            }
        }
        if (N <= m_MaxPrime) {
            const uint64_t val = N.get_ui();
            return m_SmallPrimes[val >> 1];
        }
        return false;
    }

    template <typename T>
    const bool
    Check(
        const T& N
    ) const {
        if (CheckSmall(N)) {
            return true;
        } else {
            if constexpr (std::is_same_v<T, mpz_class>) {
                return mpz_probab_prime_p(N.get_mpz_t(), 25) != 0;
            } else {
                mpz_class mpz_n(N);
                return mpz_probab_prime_p(mpz_n.get_mpz_t(), 25) != 0;
            }
        }
    }

    const uint64_t
    Max(void) const {
        return m_MaxPrime;
    }

private:

    void Build(
        std::span<const uint8_t> PrimeGaps
    ) {
        // Scan through to find the largest prime
        m_MaxPrime = 2;
        size_t gap_index = 1; // Start after the first gap which is for prime 2
        while (gap_index < PrimeGaps.size()) {
            // Get the next VLE-encoded gap
            uint64_t gap = 0;
            uint8_t shift = 0;
            uint8_t byte;
            do {
                if (gap_index >= PrimeGaps.size()) {
                    throw std::runtime_error("Ran out of prime gaps");
                }
                byte = PrimeGaps[gap_index];
                gap |= static_cast<uint64_t>(byte & 0x7F) << shift;
                shift += 7;
                gap_index++;
            } while ((byte & 0x80) != 0);
            m_MaxPrime += gap;
        }
        
        // Initialize the sieve. We don't store even numbers, so the number of entries is m_MaxPrime / 2 + 1
        m_SmallPrimes.resize(m_MaxPrime / 2 + 1, false);
        
        // Now mark odd primes using the gaps (skip prime 2 since it's handled separately in check())
        uint64_t prime = 2;
        gap_index = 1; // Reset to start after the first gap
        while (gap_index < PrimeGaps.size()) {
            // Get the next VLE-encoded gap
            uint64_t gap = 0;
            uint8_t shift = 0;
            uint8_t byte;
            do {
                byte = PrimeGaps[gap_index];
                gap |= static_cast<uint64_t>(byte & 0x7F) << shift;
                shift += 7;
                gap_index++;
            } while ((byte & 0x80) != 0);
            prime += gap;
            m_SmallPrimes[prime >> 1] = true;
        }
    }

    std::vector<bool> m_SmallPrimes;
    uint64_t m_MaxPrime;
};