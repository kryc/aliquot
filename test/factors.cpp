#include <gmpxx.h>
#include <gtest/gtest.h>

#include "factors.hpp"

TEST(Factors, Composite)
{
    PrimeFactors factors;
    factors.add_factor(2);
    factors.add_factor(2);
    factors.add_factor(3);
    factors.add_factor(5);
    auto composites = factors.get_composite(true);
    // The composites of 2^2 * 3^1 * 5^1 are:
    // 1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60
    std::vector<mpz_class> expected = {
        1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60
    };
    ASSERT_EQ(composites.size(), expected.size());
    for (size_t i = 0; i < expected.size(); ++i) {
        EXPECT_EQ(composites[i], expected[i]);
    }
}

TEST(Factors, CompositeManyFactors)
{
    PrimeFactors factors;
    // Building 2^3 * 3^2 * 5 * 7 = 2520
    factors.add_factor(2);
    factors.add_factor(2);
    factors.add_factor(2);
    factors.add_factor(3);
    factors.add_factor(3);
    factors.add_factor(5);
    factors.add_factor(7);
    auto composites = factors.get_composite(true);
    
    // 2520 has (3+1)(2+1)(1+1)(1+1) = 48 divisors
    std::vector<mpz_class> expected = {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
        12, 14, 15, 18, 20, 21, 24, 28, 30, 35,
        36, 40, 42, 45, 56, 60, 63, 70, 72, 84,
        90, 105, 120, 126, 140, 168, 180, 210, 252, 280,
        315, 360, 420, 504, 630, 840, 1260, 2520
    };
    ASSERT_EQ(composites.size(), expected.size());
    for (size_t i = 0; i < expected.size(); ++i) {
        EXPECT_EQ(composites[i], expected[i]);
    }
}

TEST(Factors, CompositeLargePrimes)
{
    PrimeFactors factors;
    // Using large primes > 2^20 (1048576)
    // 1048583^2 * 1048589 (first two primes after 2^20)
    mpz_class p1(1048583);  // First prime after 2^20
    mpz_class p2(1048589);  // Second prime after 2^20
    
    factors.add_factor(p1);
    factors.add_factor(p1);
    factors.add_factor(p2);
    auto composites = factors.get_composite(true);
    
    // Should have (2+1)(1+1) = 6 divisors
    std::vector<mpz_class> expected = {
        mpz_class("1"),
        p1,                    // 1048583
        p2,                    // 1048589
        p1 * p1,              // 1099526307889
        p1 * p2,              // 1099518288787
        p1 * p1 * p2          // 1152956817116487
    };
    ASSERT_EQ(composites.size(), expected.size());
    for (size_t i = 0; i < expected.size(); ++i) {
        EXPECT_EQ(composites[i], expected[i]);
    }
}