#include <iostream>

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
        if (sum == 0 || sum == current) {
            break;
        }
        sequence.push_back(sum);
        if (verbose)
            std::cout << sum << std::endl;
        current = sum;
    }
    return sequence;
}