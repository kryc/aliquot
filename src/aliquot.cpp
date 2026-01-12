#include <iostream>
#include <span>

#include "gmpxx.h"

#include "aliquot.hpp"
#include "primefactors.hpp"

mpz_class
sum_of_divisors(
    const mpz_class& n
)
{
    // Get prime factors of n
    auto factors = prime_factors(n);
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
    const bool verbose
)
{
    std::vector<mpz_class> sequence;
    mpz_class current = n;
    while (true) {
        mpz_class sum = sum_of_divisors(current);
        if (sum == 0) {
            break;
        }
        if (verbose)
            std::cout << sum << std::endl;
        sequence.push_back(sum);
        if (sum == current || detect_loop(sequence, sum)) {
            break;
        }
        current = sum;
    }
    return sequence;
}