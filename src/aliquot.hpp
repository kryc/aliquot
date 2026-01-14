#pragma once

#include <cinttypes>
#include <string_view>
#include <thread>
#include <tuple>
#include <vector>

#include <gmpxx.h>

#include "isprime.hpp"
#include "primefactorcache.hpp"

std::tuple<mpz_class, PrimeFactors>
SumOfDivisors(
    const mpz_class& N,
    PrimeFactorCache<>& Cache,
    const size_t NumThreads
);

std::tuple<mpz_class, PrimeFactors>
SumOfDivisors(
    const mpz_class& N
);

std::vector<mpz_class>
AliquotSequence(
    const mpz_class& N,
    const std::string_view CachePath = "",
    const bool Verbose = false,
    const size_t NumThreads = std::thread::hardware_concurrency()
);