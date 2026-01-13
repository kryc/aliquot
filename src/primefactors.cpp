#include <future>
#include <iostream>
#include <mutex>
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
    PrimeFactorCache& cache
)
{
    IsPrime is_prime;
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

bool
prime_factors_in_range(
    const mpz_class& n,
    const mpz_class& min_factor,
    const mpz_class& max_factor,
    const size_t modulus,
    std::span<const uint64_t> wheel_gaps,
    PrimeFactors& prime_factors,
    std::mutex& factor_mutex,
    bool& found
)
{
    // Take a writeable copy of n
    mpz_class target = n;

    mpz_class candidate = min_factor;
    // Ensure candidate is congruent to 1 mod modulus
    // The min_factor should already be aligned to the wheel, but just in case
    if (candidate % modulus != 0) {
        std::lock_guard<std::mutex> lock(factor_mutex);
        std::cerr << "min_factor: " << candidate << " is not aligned to wheel modulus: " << modulus << std::endl;
        throw std::runtime_error("min_factor is not aligned to wheel modulus.");
    }
    candidate += 1;

    while (candidate < max_factor && !found) {
        for (auto gapword : wheel_gaps) {
            for (size_t i = 0; i < kGapsPerWord; ++i) {
                // Check if candidate divides n
                if (candidate != 1 && mpz_divisible_p(target.get_mpz_t(), candidate.get_mpz_t())) {
                    // Lock and add factor
                    std::lock_guard<std::mutex> lock(factor_mutex);
                    while (mpz_divisible_p(target.get_mpz_t(), candidate.get_mpz_t())) {
                        prime_factors.add_factor(candidate);
                        target /= candidate;
                        // std::cout << "Found factor: " << candidate << ", Remaining target: " << target << std::endl;
                    }
                    if (target == 1) {
                        found = true;
                        return true;
                    }
                    // Check if we have found all factors
                    mpz_class current_product = prime_factors.product();
                    if (current_product >= n) {
                        found = true;
                        return true;
                    }
                    // mpz_class remainder = n / current_product;
                    // if (mpz_probab_prime_p(remainder.get_mpz_t(), 25) != 0) {
                    //     prime_factors.add_factor(remainder);
                    //     found = true;
                    //     return true;
                    // }
                }
                // Get next candidate
                const uint64_t increment = gapword & kGapMask;
                gapword >>= kBitsPerWheelGap;
                candidate += increment;
            }
        }
    }
    return false;
}

PrimeFactors
prime_factors_mt(
    const mpz_class& n,
    PrimeFactorCache& cache,
    const size_t num_threads
)
{
    // First we need to figure out which wheel modulus to use
    // It is always more efficient to use all compute cores than
    // a large modulus. So we pick the smallest modulus that divides
    // the range evenly among threads.
    mpz_class modulus = n / num_threads;
    // Round down to nearest wheel modulus
    if (modulus >= 223092870) {
        modulus = 223092870;
    } else if (modulus >= 9699690) {
        modulus = 9699690;
    } else if (modulus >= 510510) {
        modulus = 510510;
    } else if (modulus >= 30030) {
        modulus = 30030;
    } else if (modulus >= 2310) {
        modulus = 2310;
    } else if (modulus >= 210) {
        modulus = 210;
    } else {
        throw std::runtime_error("Number too small for multi-threaded factorization.");
    }

    // std::cout << "Using wheel modulus: " << modulus << std::endl;

    // Get the wheel gaps
    const size_t modulus_ui = modulus.get_ui();
    std::span<const uint64_t> wheel_gaps = get_wheel(modulus_ui);

    // Now work out the range that we want to search. We will use 0 - sqrt(n)
    mpz_class sqrt_n;
    mpz_sqrt(sqrt_n.get_mpz_t(), n.get_mpz_t());
    mpz_class min_factor = 0;
    // Round up sqrt_n to nearest multiple of modulus
    mpz_class max_factor = (sqrt_n + modulus - 1) / modulus * modulus;
    mpz_class range = max_factor - min_factor;
    mpz_class blocks = range / modulus;
    mpz_class blocksize = blocks / num_threads * modulus;
    // std::cout << "Factoring range: 0 to " << max_factor << " using " << num_threads << " threads." << std::endl;

    // Divide the search space among threads
    std::vector<std::future<bool>> futures;
    std::mutex factor_mutex;
    bool found = false;

    PrimeFactors local_factors;
    // Divide out the small primes below the wheel modulus
    mpz_class remainder = n;
    for (auto& prime : get_primes_for_wheel_modulus(modulus.get_ui())) {
        while (remainder % prime == 0) {
            local_factors.add_factor(prime);
            remainder /= prime;
        }
    }
    
    for (size_t i = 0; i < num_threads; ++i) {
        mpz_class thread_min = min_factor + i * blocksize;
        mpz_class thread_max = (i == num_threads - 1) ? max_factor : (thread_min + blocksize - 1);

        futures.push_back(std::async(std::launch::async, [t_min = thread_min, t_max = thread_max, &remainder, modulus_ui, wheel_gaps, &local_factors, &factor_mutex, &found]() {
            return prime_factors_in_range(remainder, t_min, t_max, modulus_ui, wheel_gaps, local_factors, factor_mutex, found);
        }));
    }

    // PrimeFactors total_factors;
    for (auto& fut : futures) {
        fut.get();
    }

    const mpz_class product = local_factors.product();
    if (product == n) {
        return local_factors;
    }
    else if (product > n) {
        throw std::runtime_error("Product of found factors exceeds n.");
    }
    else if (product < n) {
        // If we didn't find all factors, check if the remaining quotient is prime
        const mpz_class rem = n / product;
        if (rem > 1 && mpz_probab_prime_p(rem.get_mpz_t(), 25) != 0) {
                local_factors.add_factor(rem);
        } else {
            throw std::runtime_error("Failed to fully factor the number in the given range.");
        }
    }

    return local_factors;
}

PrimeFactors
prime_factors(
    const mpz_class& n,
    PrimeFactorCache& cache,
    const size_t num_threads
)
{
     // If the number is small, use the linear method
    if (n < 1'000'000) {
        return prime_factors_linear(n, cache);
    }

    // For larger numbers, use multi-threaded factorization
    return prime_factors_mt(n, cache, num_threads);
}

PrimeFactors
prime_factors(
    const mpz_class& n
)
{
    PrimeFactorCache cache("");
    const size_t num_threads = std::thread::hardware_concurrency();
    return prime_factors(n, cache, num_threads);
}