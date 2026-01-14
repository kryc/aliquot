#pragma once

#include <thread>
#include <vector>

#include <gmpxx.h>

#include "primefactorcache.hpp"
#include "factors.hpp"
#include "isprime.hpp"

PrimeFactors
PrimeFactorsLinear(
    const mpz_class& N,
    PrimeFactorCache<>& Cache
);

PrimeFactors
PrimeFactorsInRange(
    const mpz_class& N,
    const mpz_class& MinFactor,
    const mpz_class& MaxFactor
);

PrimeFactors
PrimeFactorsMT(
    const mpz_class& N,
    PrimeFactorCache<>& Cache,
    const size_t NumThreads = std::thread::hardware_concurrency()
);

PrimeFactors
GetPrimeFactors(
    const mpz_class& N,
    PrimeFactorCache<>& Cache,
    const size_t NumThreads = std::thread::hardware_concurrency()
);

PrimeFactors
GetPrimeFactors(
    const mpz_class& N
);