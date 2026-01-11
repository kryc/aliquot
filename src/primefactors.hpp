#pragma once

#include <vector>

#include <gmpxx.h>

#include "factors.hpp"

PrimeFactors
prime_factors(
    const mpz_class& n
);