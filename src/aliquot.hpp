#pragma once

#include <cinttypes>
#include <string_view>
#include <vector>

#include <gmpxx.h>

mpz_class
sum_of_divisors(
    const mpz_class& n
);

std::vector<mpz_class>
aliquot_sequence(
    const mpz_class& n,
    const std::string_view cache_path = "",
    const bool verbose = false
);