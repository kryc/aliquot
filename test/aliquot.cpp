#include <gtest/gtest.h>

#include <gmpxx.h>

#include "aliquot.hpp"
#include "primefactors.hpp"

TEST(Aliquot, SumOfDivisors)
{
    mpz_class n = 10; // 1, 2, 5
    auto sum = sum_of_divisors(n);
    EXPECT_EQ(sum, 8); // 1 + 2 + 5 = 8
    n = 8;
    sum = sum_of_divisors(n); // 1, 2, 4
    EXPECT_EQ(sum, 7); // 1 + 2 + 4 = 7
}

TEST(Aliquot, AliquotSequence)
{
    mpz_class n = 12; // 1, 2, 3, 4, 6
    auto sequence = aliquot_sequence(n);
    std::vector<mpz_class> expected = {16, 15, 9, 4, 3, 1};
    ASSERT_EQ(sequence.size(), expected.size());
    for (size_t i = 0; i < expected.size(); ++i) {
        EXPECT_EQ(sequence[i], expected[i]);
    }
}