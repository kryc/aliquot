#pragma once

#include <thread>
#include <vector>

#include <gmpxx.h>

#include "primefactorcache.hpp"
#include "factors.hpp"
#include "isprime.hpp"

PrimeFactors
prime_factors_linear(
    const mpz_class& n,
    PrimeFactorCache& cache
);

PrimeFactors
prime_factors_in_range(
    const mpz_class& n,
    const mpz_class& min_factor,
    const mpz_class& max_factor
);

PrimeFactors
prime_factors_mt(
    const mpz_class& n,
    PrimeFactorCache& cache,
    const size_t num_threads = std::thread::hardware_concurrency()
);

PrimeFactors
prime_factors(
    const mpz_class& n,
    PrimeFactorCache& cache,
    const size_t num_threads = std::thread::hardware_concurrency()
);

PrimeFactors
prime_factors(
    const mpz_class& n
);