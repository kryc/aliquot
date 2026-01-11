#pragma once

#include <algorithm>
#include <cinttypes>
#include <map>
#include <vector>

#include <gmpxx.h>

// A placeholder for an efficient prime factor storage structure
class PrimeFactors {
public:
    void add_factor(
        const mpz_class& factor
    ) {
        factor_counts[factor]++;
    }

    bool has_factor(
        const mpz_class& factor
    ) const {
        return factor_counts.find(factor) != factor_counts.end();
    }

    size_t count(
        void
    ) const {
        size_t total = 0;
        for (const auto& pair : factor_counts) {
            total += pair.second;
        }
        return total;
    }

    size_t size(
        void
    ) const {
        return factor_counts.size();
    }

    size_t count_of(
        const mpz_class& factor
    ) const {
        auto it = factor_counts.find(factor);
        if (it != factor_counts.end()) {
            return it->second;
        }
        return 0;
    }

    bool empty(
        void
    ) const {
        return factor_counts.empty();
    }

    std::vector<std::pair<mpz_class, size_t>> to_vector() const {
        std::vector<std::pair<mpz_class, size_t>> vec;
        for (const auto& pair : factor_counts) {
            vec.emplace_back(pair.first, pair.second);
        }
        return vec;
    }

    // Function to convert prime factors to a vector of composite factors
    // By iterating over all combinations of prime factor powers
    std::vector<mpz_class>
    get_composite(
        const bool sorted = false
    ) const {
        std::vector<mpz_class> composites;
        composites.push_back(1); // Start with 1 as the first composite
        for (const auto& [prime, count] : factor_counts) {
            size_t current_size = composites.size();
            mpz_class prime_power = 1;
            for (size_t i = 1; i <= count; ++i) {
                prime_power *= prime;
                for (size_t j = 0; j < current_size; ++j) {
                    composites.push_back(composites[j] * prime_power);
                }
            }
        }
        if (sorted)
        {
            std::sort(composites.begin(), composites.end());
        }
        return composites;
    }

    std::vector<uint8_t>
    serialize(
        void
    ) const {
        std::vector<uint8_t> data;
        return data;
    }

    mpz_class
    max_factor(
        void
    ) const {
        if (factor_counts.empty()) {
            return 0;
        }
        return std::max_element(
            factor_counts.begin(),
            factor_counts.end(),
            [](const auto& a, const auto& b) {
                return a.first < b.first;
            })->first;
    }

    size_t
    max_power(
        void
    ) const {
        if (factor_counts.empty()) {
            return 0;
        }
        return std::max_element(
            factor_counts.begin(),
            factor_counts.end(),
            [](const auto& a, const auto& b) {
                return a.second < b.second;
            })->second;
    }
private:
    std::map<mpz_class, size_t> factor_counts;
};