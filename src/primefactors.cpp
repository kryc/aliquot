#include <future>
#include <iostream>
#include <thread>
#include <vector>

#include <gmpxx.h>

#include "factors.hpp"
#include "isprime.hpp"
#include "primefactorcache.hpp"
#include "primefactors.hpp"
#include "primes.hpp"

PrimeFactors
prime_factors_linear(
    const mpz_class& n,
    PrimeFactorCache& cache,
    IsPrime& is_prime
)
{
    // Get the span of prime gaps
    auto gaps = get_prime_gaps();

    PrimeFactors prime_factors;

    std::vector<mpz_class> factors;
    mpz_class prime = 2;
    size_t gap_index = 1;
    mpz_class remainder = n;

    while (remainder > 1 && gap_index < gaps.size()) {
        // Check if the remainder is prime
        if (is_prime.check_small(remainder)) {
            prime_factors.add_factor(remainder);
            return prime_factors;
        }
        // Check if we have the factor in cache
        // auto cached_factors = cache.product_exists(remainder.get_ui());
        // if (cached_factors.has_value()) {
        //     auto cached = cached_factors.value();
        //     prime_factors.update(cached);
        //     return prime_factors;
        // }
        // Check if prime divides remainder
        while (remainder % prime == 0) {
            prime_factors.add_factor(prime);
            remainder /= prime;
        }

        // Get the next VLE-encoded gap
        uint64_t gap = 0;
        uint8_t shift = 0;
        uint8_t byte;
        do {
            byte = gaps[gap_index];
            gap |= static_cast<uint64_t>(byte & 0x7F) << shift;
            shift += 7;
            gap_index++;
        } while ((byte & 0x80) != 0);
        prime += gap;
    }

    if (remainder != 1) {
        // std::cout << "Continuing with prime " << prime << " for remainder " << remainder << std::endl;
    }

    // If remainder is not 1, then we need to continue
    while (remainder > 1) {
        // Check if prime divides remainder
        while (remainder % prime == 0) {
            prime_factors.add_factor(prime);
            remainder /= prime;
        }
        // Get next prime
        mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
    }

    return prime_factors;
}

PrimeFactors
prime_factors_in_range(
    const mpz_class& n,
    const mpz_class& min_factor,
    const mpz_class& max_factor,
    IsPrime& is_prime
)
{
    PrimeFactors prime_factors;

    mpz_class factor = min_factor;
    if (!is_prime.check(factor)) {
        mpz_nextprime(factor.get_mpz_t(), factor.get_mpz_t());
    }

    while (factor <= max_factor) {
        // Check if factor divides n (not a local remainder)
        mpz_class temp = n;
        while (temp % factor == 0) {
            prime_factors.add_factor(factor);
            temp /= factor;
        }
        // Get next prime factor
        mpz_nextprime(factor.get_mpz_t(), factor.get_mpz_t());
    }

    return prime_factors;
}

PrimeFactors
prime_factors_mt(
    const mpz_class& n,
    PrimeFactorCache& cache,
    IsPrime& is_prime,
    const size_t num_threads
)
{
    // The minimum prime factor we need to check is 2
    // The maximum prime factor we need to check is n/2
    mpz_class min_factor = 2;
    mpz_class max_factor = n / 2;

    // Divide the search space among threads
    std::vector<std::future<PrimeFactors>> futures;
    mpz_class range_size = (max_factor - min_factor) / num_threads;

    for (size_t i = 0; i < num_threads; ++i) {
        mpz_class thread_min = min_factor + i * range_size;
        mpz_class thread_max = (i == num_threads - 1) ? max_factor : (thread_min + range_size - 1);
        // std::cout << "Thread " << i << " searching range [" << thread_min << ", " << thread_max << "]" << std::endl;

        futures.push_back(std::async(std::launch::async, [thread_min, thread_max, &n, &is_prime]() {
            return prime_factors_in_range(n, thread_min, thread_max, is_prime);
        }));
    }

    PrimeFactors total_factors;
    for (auto& fut : futures) {
        PrimeFactors thread_factors = fut.get();
        total_factors.update(thread_factors);
    }

    if (total_factors.product() < n) {
        throw std::runtime_error("Failed to fully factor the number in the given range");
    }

    return total_factors;
}

PrimeFactors
prime_factors(
    const mpz_class& n,
    PrimeFactorCache& cache,
    IsPrime& is_prime,
    const size_t num_threads
)
{
     // If the number is small, use the linear method
    if (n < 1'000'000) {
        return prime_factors_linear(n, cache, is_prime);
    }

    // For larger numbers, use multi-threaded factorization
    return prime_factors_mt(n, cache, is_prime, num_threads);
}