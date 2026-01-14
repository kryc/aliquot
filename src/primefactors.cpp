#include <atomic>
#include <bit>
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
PrimeFactorsLinear(
    const mpz_class& N,
    PrimeFactorCache<>& Cache
)
{
    IsPrime prime_checker;
    // Get the span of prime gaps
    auto gaps = GetPrimeGaps();

    PrimeFactors prime_factors;

    std::vector<mpz_class> factors;
    uint64_t prime = 2;
    size_t gap_index = 1;
    mpz_class remainder = N;

    while (remainder > 1 && gap_index < gaps.size()) {
        // Check if the remainder is prime
        if (prime_checker.CheckSmall(remainder)) {
            prime_factors.AddFactor(remainder);
            return prime_factors;
        }
        // Check if we have the factor in cache
        // auto cached_factors = cache.ProductExists(remainder);
        // if (cached_factors.has_value()) {
        //     auto cached = cached_factors.value();
        //     prime_factors.Update(cached);
        //     return prime_factors;
        // }
        // Check if prime divides remainder
        while (mpz_divisible_ui_p(remainder.get_mpz_t(), prime)) {
            prime_factors.AddFactor(prime);
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

    // If remainder is not 1, then we need to continue factoring
    // We use the 30-wheel to avoid using nextprime
    if (remainder != 1) {
        // Get a smallish wheel modulus
        uint32_t wheel = kWheel30;
        // Set our candidate to the last prime used and round down so it
        // is congruent to 1 mod modulus
        mpz_class candidate = prime;
        while (candidate % 30 != 1) {
            candidate -= 1;
        }
        // Now use the wheel to factor the remainder
        while (remainder > 1) {
            while (mpz_divisible_p(remainder.get_mpz_t(), candidate.get_mpz_t())) {
                prime_factors.AddFactor(candidate);
                remainder /= candidate;
            }
            const uint32_t increment = wheel & kWheel30Mask;
            wheel = std::rotr(wheel, kWheel30BitsPerGap);
            candidate += increment;
        }
    }

    return prime_factors;
}

bool
PrimeFactorsInRange(
    const mpz_class& N,
    const mpz_class& MinFactor,
    const mpz_class& MaxFactor,
    const IsPrime& PrimeChecker,
    const size_t Modulus,
    std::span<const uint64_t> WheelGaps,
    PrimeFactors& FoundFactors,
    std::mutex& Mutex,
    std::atomic<bool>& Found
)
{
    mpz_class candidate = MinFactor;
    // Ensure candidate is congruent to 1 mod modulus
    // The min_factor should already be aligned to the wheel, but just in case
    if (candidate % Modulus != 0) {
        std::lock_guard<std::mutex> lock(Mutex);
        throw std::runtime_error("min_factor is not aligned to wheel modulus.");
    }
    candidate += 1;

    while (candidate < MaxFactor && !Found.load()) {
        for (auto gapword : WheelGaps) {
            for (size_t i = 0; i < kGapsPerWord; ++i) {
                // Check if candidate divides n
                if (candidate != 1 && mpz_divisible_p(N.get_mpz_t(), candidate.get_mpz_t())) {
                    // Lock and add factor
                    std::lock_guard<std::mutex> lock(Mutex);
                    // Calculate the remaining quotient that hasn't been factored yet
                    mpz_class current_product = FoundFactors.Product();
                    mpz_class quotient = N / current_product;
                    // if (!mpz_divisible_p(quotient.get_mpz_t(), candidate.get_mpz_t())) {
                    //     std::cerr << "Error: Candidate " << candidate << " does not divide current quotient " << quotient << std::endl;
                    //     throw std::runtime_error("Candidate does not divide the current quotient.");
                    // }
                    // Add all powers of this factor that divide the quotient
                    while (mpz_divisible_p(quotient.get_mpz_t(), candidate.get_mpz_t())) {
                        FoundFactors.AddFactor(candidate);
                        quotient /= candidate;
                    }

                    if (quotient == N) {
                        Found.store(true);
                        return true;
                    } else if (quotient > N) {
                        throw std::runtime_error("Product of found factors exceeds n.");
                    } else {
                        // See if we can return early if the remaining quotient is prime
                        mpz_class remainder = N / quotient;
                        if (PrimeChecker.Check(remainder)) {
                            FoundFactors.AddFactor(remainder);
                            Found.store(true);
                            return true;
                        }
                    }
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
PrimeFactorsMT(
    const mpz_class& N,
    PrimeFactorCache<>& Cache,
    const size_t NumThreads
)
{
    IsPrime prime_checker;
    // First we work out the range that we want to search. We will use 0 - sqrt(n)
    mpz_class sqrt_n;
    mpz_sqrt(sqrt_n.get_mpz_t(), N.get_mpz_t());

    // Next we need to figure out which wheel modulus to use
    // It is always more efficient to use all compute cores than
    // a large modulus. So we pick the smallest modulus that divides
    // the range evenly among threads.
    mpz_class modulus = sqrt_n / NumThreads;
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
    } else if (modulus >= 30) {
        modulus = 30;
    } else {
        std::cerr << "Modulus: " << modulus << " is too small for wheel factorization." << std::endl;
        throw std::runtime_error("Number too small for multi-threaded factorization.");
    }

    // Get the wheel gaps
    const size_t modulus_ui = modulus.get_ui();
    std::span<const uint64_t> wheel_gaps = GetWheel(modulus_ui);
    
    // Round up sqrt_n to nearest multiple of modulus
    mpz_class max_factor = (sqrt_n + modulus - 1) / modulus * modulus;
    // std::cout << "Factoring range: 0 to " << max_factor << " using " << NumThreads << " threads." << std::endl;

    // Divide the search space among threads using interleaved distribution
    // Each thread processes every Nth block (where N = NumThreads)
    std::vector<std::future<bool>> futures;
    std::mutex factor_mutex;
    std::atomic<bool> found = false;

    PrimeFactors local_factors;
    // Divide out the small primes below the wheel modulus
    mpz_class remainder = N;
    for (auto& prime : GetPrimesForWheelModulus(modulus.get_ui())) {
        while (remainder % prime == 0) {
            local_factors.AddFactor(prime);
            remainder /= prime;
        }
    }
    
    // Launch threads with interleaved block distribution
    for (size_t i = 0; i < NumThreads; ++i) {
        futures.push_back(std::async(std::launch::async, [thread_id = i, NumThreads, modulus, &max_factor, &N, modulus_ui, wheel_gaps, &local_factors, &factor_mutex, &found, &prime_checker]() {
            // This thread processes blocks: thread_id, thread_id + num_threads, thread_id + 2*num_threads, ...
            mpz_class block_start = thread_id * modulus;
            while (block_start < max_factor && !found) {
                mpz_class block_end = block_start + modulus;
                if (block_end > max_factor) {
                    block_end = max_factor;
                }
                if (PrimeFactorsInRange(N, block_start, block_end, prime_checker, modulus_ui, wheel_gaps, local_factors, factor_mutex, found)) {
                    return true;
                }
                block_start += NumThreads * modulus;
            }
            return false;
        }));
    }

    // PrimeFactors total_factors;
    for (auto& fut : futures) {
        fut.get();
    }

    const mpz_class product = local_factors.Product();
    if (product == N) {
        return local_factors;
    } else if (product > N) {
        throw std::runtime_error("Product of found factors exceeds n.");
    } else if (product < N) {
        // If we didn't find all factors, check if the remaining quotient is prime
        const mpz_class rem = N / product;
        if (rem > 1 && mpz_probab_prime_p(rem.get_mpz_t(), 25) != 0) {
                local_factors.AddFactor(rem);
        } else {
            throw std::runtime_error("Failed to fully factor the number in the given range.");
        }
    }

    return local_factors;
}

PrimeFactors
GetPrimeFactors(
    const mpz_class& N,
    PrimeFactorCache<>& Cache,
    const size_t NumThreads
)
{
    auto cached = Cache.ProductExists(N.get_ui());
    if (cached.has_value()) {
        return cached.value();
    }
     // If the number is small, use the linear method
    if (N < 3'000'000) {
        return PrimeFactorsLinear(N, Cache);
    }

    // For larger numbers, use multi-threaded factorization
    return PrimeFactorsMT(N, Cache, NumThreads);
}

PrimeFactors
GetPrimeFactors(
    const mpz_class& N
)
{
    PrimeFactorCache<> cache("");
    const size_t num_threads = std::thread::hardware_concurrency();
    return GetPrimeFactors(N, cache, num_threads);
}