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
        build_(get_prime_gaps());
    }
    IsPrime(std::span<const uint8_t> prime_gaps) {
        build_(prime_gaps);
    }

    template <typename T>
    const bool
    check_small(
        const T& n
    ) const {
        if (n < 2) {
            return false;
        }
        if (n == 2) {
            return true;
        }
        if constexpr (std::is_same_v<T, mpz_class>) {
            if ((n.get_ui() & 1) == 0) {
                return false;  // Even numbers > 2 are not prime
            }
        } else {
            if ((n & 1) == 0) {
                return false;  // Even numbers > 2 are not prime
            }
        }
        if (n <= max_prime_) {
            const uint64_t val = n.get_ui();
            return small_primes_[val >> 1];
        }
        return false;
    }

    template <typename T>
    const bool
    check(
        const T& n
    ) const {
        if (check_small(n)) {
            return true;
        } else {
            if constexpr (std::is_same_v<T, mpz_class>) {
                return mpz_probab_prime_p(n.get_mpz_t(), 25) != 0;
            } else {
                mpz_class mpz_n(n);
                return mpz_probab_prime_p(mpz_n.get_mpz_t(), 25) != 0;
            }
        }
    }

    const uint64_t
    max(void) const {
        return max_prime_;
    }

private:

    void build_(
        std::span<const uint8_t> prime_gaps
    ) {
        // Scan through to find the largest prime
        max_prime_ = 2;
        size_t gap_index = 1; // Start after the first gap which is for prime 2
        while (gap_index < prime_gaps.size()) {
            // Get the next VLE-encoded gap
            uint64_t gap = 0;
            uint8_t shift = 0;
            uint8_t byte;
            do {
                if (gap_index >= prime_gaps.size()) {
                    throw std::runtime_error("Ran out of prime gaps");
                }
                byte = prime_gaps[gap_index];
                gap |= static_cast<uint64_t>(byte & 0x7F) << shift;
                shift += 7;
                gap_index++;
            } while ((byte & 0x80) != 0);
            max_prime_ += gap;
        }
        
        // Initialize the sieve. We don't store even numbers, so the number of entries is max_prime_ / 2 + 1
        small_primes_.resize(max_prime_ / 2 + 1, false);
        
        // Now mark odd primes using the gaps (skip prime 2 since it's handled separately in check())
        uint64_t prime = 2;
        gap_index = 1; // Reset to start after the first gap
        while (gap_index < prime_gaps.size()) {
            // Get the next VLE-encoded gap
            uint64_t gap = 0;
            uint8_t shift = 0;
            uint8_t byte;
            do {
                byte = prime_gaps[gap_index];
                gap |= static_cast<uint64_t>(byte & 0x7F) << shift;
                shift += 7;
                gap_index++;
            } while ((byte & 0x80) != 0);
            prime += gap;
            small_primes_[prime >> 1] = true;
        }
    }

    std::vector<bool> small_primes_;
    uint64_t max_prime_;
};