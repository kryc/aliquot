#pragma once

#include <vector>

#include <gmpxx.h>

mpz_class
sum_of_divisors(
    const mpz_class& n
);

std::vector<mpz_class>
aliquot_sequence(
    const mpz_class& n,
    const bool verbose = false
);