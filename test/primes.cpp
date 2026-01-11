#include <gtest/gtest.h>

#include "primes.hpp"
#include "primefactors.hpp"

TEST(Primes, GeneratePrimes) {
    auto gaps = generate_prime_gaps(100, false);
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
    auto gaps_opt = get_prime_gaps(256ull);
    if (gaps_opt.has_value()) {
        auto gaps = gaps_opt.value();
        ASSERT_FALSE(gaps.empty());
        // Just check the first few gaps
        EXPECT_TRUE(gaps[0] == 2);  // gap 2 (0 to 2)
        EXPECT_TRUE(gaps[1] == 1);  // gap 1 (2 to 3)
        EXPECT_TRUE(gaps[2] == 2);  // gap 2 (3 to 5)
    } else {
        GTEST_SKIP() << "No prime gaps loaded.";
    }
}

TEST(Primes, PrimeFactorsSmall)
{
    mpz_class n = 100;
    auto factors = prime_factors(n);
    EXPECT_EQ(factors.size(), 2); // 2^2 * 5^2
    EXPECT_EQ(factors.count(), 4);
    EXPECT_EQ(factors.count_of(2), 2);
    EXPECT_EQ(factors.count_of(5), 2);
}

TEST(Primes, PrimeFactorsMedium)
{
    mpz_class n(65536); // 2^16
    auto factors = prime_factors(n);
    EXPECT_EQ(factors.size(), 1);
    EXPECT_EQ(factors.count(), 16);
    EXPECT_EQ(factors.count_of(2), 16);
}

// Test the overflow logic in prime factorization
// where we run out of precomputed prime gaps
// We test with a multiple of the first prime after 2^16, which is 65537
TEST(Primes, PrimeFactorsLarge)
{
    mpz_class n(131074); // 65537^1 * 2^1
    auto factors = prime_factors(n);
    EXPECT_EQ(factors.size(), 2);
    EXPECT_EQ(factors.count(), 2);
    EXPECT_EQ(factors.count_of(2), 1);
    EXPECT_EQ(factors.count_of(65537), 1);
}

TEST(Primes, GetNthPrime)
{
    // Test first 10 primes
    std::vector<mpz_class> expected_primes = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29
    };
    for (size_t i = 0; i < expected_primes.size(); ++i) {
        auto prime = get_nth_prime(i);
        EXPECT_EQ(prime, expected_primes[i]);
    }
    // Test with a very large (bigger than precomputed gaps) index
    size_t large_index = 70000; // 70000th prime is 882377
    auto large_prime = get_nth_prime(large_index);
    EXPECT_EQ(large_prime, mpz_class(882389));
}

TEST(Primes, GetPrimeIndex)
{
    // Test first 10 primes
    std::vector<mpz_class> primes = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29
    };
    for (size_t i = 0; i < primes.size(); ++i) {
        auto index = get_prime_index(primes[i]);
        EXPECT_EQ(index, i);
    }
    // Test with a very large prime (bigger than precomputed gaps)
    mpz_class large_prime(882389); // 70000th prime
    auto large_index = get_prime_index(large_prime);
    EXPECT_EQ(large_index, 70000);
}