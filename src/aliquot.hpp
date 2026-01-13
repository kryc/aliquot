#pragma once

#include <cinttypes>
#include <string_view>
#include <thread>
#include <vector>

#include <gmpxx.h>

#include "isprime.hpp"
#include "primefactorcache.hpp"

mpz_class
sum_of_divisors(
    const mpz_class& n,
    PrimeFactorCache& cache,
    IsPrime& is_prime,
    const size_t num_threads
);

std::vector<mpz_class>
aliquot_sequence(
    const mpz_class& n,
    const std::string_view cache_path = "",
    const bool verbose = false,
    const size_t num_threads = std::thread::hardware_concurrency()
);