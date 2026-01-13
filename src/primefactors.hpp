#pragma once

#include <vector>

#include <gmpxx.h>

#include "primefactorcache.hpp"
#include "factors.hpp"
#include "isprime.hpp"

PrimeFactors
prime_factors(
    const mpz_class& n,
    PrimeFactorCache& cache,
    IsPrime& is_prime
);