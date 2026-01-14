#include <gtest/gtest.h>

#include "isprime.hpp"
#include "primes.hpp"
#include "primefactors.hpp"

TEST(Primes, GeneratePrimes) {
    auto gaps = GeneratePrimeGaps(100, false);
    ASSERT_FALSE(gaps.empty());
    // The first few primes are 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97
    EXPECT_TRUE(gaps[0] == 2);  // gap 2 (0 to 2)
    EXPECT_TRUE(gaps[1] == 1);  // gap 1 (2 to 3)
    EXPECT_TRUE(gaps[2] == 2);  // gap 2 (3 to 5)
    EXPECT_TRUE(gaps[3] == 2);  // gap 2 (5 to 7)
    EXPECT_TRUE(gaps[4] == 4);  // gap 4 (7 to 11)
    EXPECT_TRUE(gaps[5] == 2);  // gap 2 (11 to 13)
}

TEST(Primes, GetPrimeGaps)
{
    // This will trigger a fallback, so use a small limit
    auto gaps = GetPrimeGaps(256ull);
    ASSERT_FALSE(gaps.empty());
    // Just check the first few gaps
    EXPECT_TRUE(gaps[0] == 2);  // gap 2 (0 to 2)
    EXPECT_TRUE(gaps[1] == 1);  // gap 1 (2 to 3)
    EXPECT_TRUE(gaps[2] == 2);  // gap 2 (3 to 5)
}

TEST(Primes, PrimeFactorsSmall)
{
    mpz_class n = 100;
    auto factors = GetPrimeFactors(n);
    EXPECT_EQ(factors.Size(), 2); // 2^2 * 5^2
    EXPECT_EQ(factors.count(), 4);
    EXPECT_EQ(factors.CountOf(2), 2);
    EXPECT_EQ(factors.CountOf(5), 2);
}

TEST(Primes, PrimeFactorsMedium)
{
    mpz_class n(65536); // 2^16
    auto factors = GetPrimeFactors(n);
    EXPECT_EQ(factors.Size(), 1);
    EXPECT_EQ(factors.count(), 16);
    EXPECT_EQ(factors.CountOf(2), 16);
}

// Test the overflow logic in prime factorization
// where we run out of precomputed prime gaps
// We test with a multiple of the first prime after 2^16, which is 65537
TEST(Primes, PrimeFactorsLarge)
{
    mpz_class n(131074); // 65537^1 * 2^1
    auto factors = GetPrimeFactors(n);
    EXPECT_EQ(factors.Size(), 2);
    EXPECT_EQ(factors.count(), 2);
    EXPECT_EQ(factors.CountOf(2), 1);
    EXPECT_EQ(factors.CountOf(65537), 1);
}

// Test with a prime much larger than precomputed gaps
// close to the value 70,000. The 70000th prime is 882377
TEST(Primes, PrimeFactorsVeryLarge)
{
    mpz_class n = mpz_class(882377); // prime
    n *= 2;
    n *= 3;
    auto factors = GetPrimeFactors(n);
    EXPECT_EQ(factors.Size(), 3);
    EXPECT_EQ(factors.count(), 3);
    EXPECT_EQ(factors.CountOf(2), 1);
    EXPECT_EQ(factors.CountOf(3), 1);
    EXPECT_EQ(factors.CountOf(882377), 1);
}

TEST(Primes, GetNthPrime)
{
    // Test first 10 primes
    std::vector<mpz_class> expected_primes = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29
    };
    for (size_t i = 0; i < expected_primes.size(); ++i) {
        auto prime = GetNthPrime(i);
        EXPECT_EQ(prime, expected_primes[i]);
    }
    // Test with a very large (bigger than precomputed gaps) index
    size_t large_index = 70000; // 70000th prime is 882377
    auto large_prime = GetNthPrime(large_index);
    EXPECT_EQ(large_prime, mpz_class(882389));
}

TEST(Primes, GetPrimeIndex)
{
    // Test first 10 primes
    std::vector<mpz_class> primes = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29
    };
    for (size_t i = 0; i < primes.size(); ++i) {
        auto index = GetPrimeIndex(primes[i]);
        EXPECT_EQ(index, i);
    }
    // Test with a very large prime (bigger than precomputed gaps)
    mpz_class large_prime(882389); // 70000th prime
    auto large_index = GetPrimeIndex(large_prime);
    EXPECT_EQ(large_index, 70000);
}

TEST(IsPrime, SmallPrimes)
{
    auto gaps = GeneratePrimeGaps(100, false);
    IsPrime is_prime(gaps);
    std::vector<mpz_class> small_primes = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47
    };
    for (const auto& prime : small_primes) {
        EXPECT_TRUE(is_prime.check(prime)) << "Failed for prime: " << prime;
    }
    std::vector<mpz_class> small_non_primes = {
        0, 1, 4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20
    };
    for (const auto& non_prime : small_non_primes) {
        EXPECT_FALSE(is_prime.check(non_prime)) << "Failed for non-prime: " << non_prime;
    }
}

TEST(IsPrime, LargePrimes)
{
    auto gaps = GeneratePrimeGaps(100000, false);
    IsPrime is_prime(gaps);
    std::vector<mpz_class> large_primes = {
        999983, 1000003, 1000033, 1000037, 1000039
    };
    for (const auto& prime : large_primes) {
        EXPECT_TRUE(is_prime.check(prime)) << "Failed for prime: " << prime;
    }
    std::vector<mpz_class> large_non_primes = {
        1000000, 1000010, 1000020, 1000030
    };
    for (const auto& non_prime : large_non_primes) {
        EXPECT_FALSE(is_prime.check(non_prime)) << "Failed for non-prime: " << non_prime;
    }
}

TEST(Prime, WheelGeneration210)
{
    auto wheel = GetWheel(210);
    EXPECT_FALSE(wheel.empty());
    // The wheel should start with 1 and end at 211
    mpz_class current = 1;
    for (auto gapword : wheel) {
        for (size_t i = 0; i < kGapsPerWord; ++i) {
            const uint64_t gap = gapword & kGapMask;
            gapword >>= kBitsPerWheelGap;
            current += gap;
        }
    }
    EXPECT_EQ(current, mpz_class(211));
}

TEST(Prime, WheelGeneration2310)
{
    auto wheel = GetWheel(2310);
    EXPECT_FALSE(wheel.empty());
    // The wheel should start with 1 and end at 2311
    mpz_class current = 1;
    for (auto gapword : wheel) {
        for (size_t i = 0; i < kGapsPerWord; ++i) {
            const uint64_t gap = gapword & kGapMask;
            gapword >>= kBitsPerWheelGap;
            current += gap;
        }
    }
    EXPECT_EQ(current, mpz_class(2311));
}

TEST(Prime, WheelGeneration510510)
{
    auto wheel = GetWheel(510510);
    EXPECT_FALSE(wheel.empty());
    // The wheel should start with 1 and end at 510511
    mpz_class current = 1;
    for (auto gapword : wheel) {
        for (size_t i = 0; i < kGapsPerWord; ++i) {
            const uint64_t gap = gapword & kGapMask;
            gapword >>= kBitsPerWheelGap;
            current += gap;
        }
    }
    EXPECT_EQ(current, mpz_class(510511));
}

TEST(Prime, WheelGeneration30030)
{
    auto wheel = GetWheel(30030);
    EXPECT_FALSE(wheel.empty());
    // The wheel should start with 1 and end at 30031
    mpz_class current = 1;
    for (auto gapword : wheel) {
        for (size_t i = 0; i < kGapsPerWord; ++i) {
            const uint64_t gap = gapword & kGapMask;
            gapword >>= kBitsPerWheelGap;
            current += gap;
        }
    }
    EXPECT_EQ(current, mpz_class(30031));
}

TEST(Prime, WheelGeneration9699690)
{
    auto wheel = GetWheel(9699690);
    EXPECT_FALSE(wheel.empty());
    // The wheel should start with 1 and end at 9699691
    mpz_class current = 1;
    for (auto gapword : wheel) {
        for (size_t i = 0; i < kGapsPerWord; ++i) {
            const uint64_t gap = gapword & kGapMask;
            gapword >>= kBitsPerWheelGap;
            current += gap;
        }
    }
    EXPECT_EQ(current, mpz_class(9699691));
}

TEST(Prime, WheelGeneration223092870)
{
    auto wheel = GetWheel(223092870);
    EXPECT_FALSE(wheel.empty());
    // The wheel should start with 1 and end at 223092871
    mpz_class current = 1;
    for (auto gapword : wheel) {
        for (size_t i = 0; i < kGapsPerWord; ++i) {
            const uint64_t gap = gapword & kGapMask;
            gapword >>= kBitsPerWheelGap;
            current += gap;
        }
    }
    EXPECT_EQ(current, mpz_class(223092871));
}