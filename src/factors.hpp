#pragma once

#include <algorithm>
#include <cinttypes>
#include <map>
#include <vector>

#include <gmpxx.h>

// A placeholder for an efficient prime factor storage structure
class PrimeFactors {
public:
    void AddFactor(
        const mpz_class& Factor
    ) {
        m_FactorCounts[Factor]++;
    }

    void Update(
        const PrimeFactors& Other
    ) {
        for (const auto& [prime, count] : Other.m_FactorCounts) {
            m_FactorCounts[prime] += count;
        }
    }

    bool HasFactor(
        const mpz_class& Factor
    ) const {
        return m_FactorCounts.find(Factor) != m_FactorCounts.end();
    }

    mpz_class
    LargestFactor(
        void
    ) const {
        if (m_FactorCounts.empty()) {
            return 0;
        }
        return std::max_element(
            m_FactorCounts.begin(),
            m_FactorCounts.end(),
            [](const auto& a, const auto& b) {
                return a.first < b.first;
            })->first;
    }

    size_t Count(
        void
    ) const {
        size_t total = 0;
        for (const auto& pair : m_FactorCounts) {
            total += pair.second;
        }
        return total;
    }

    size_t Size(
        void
    ) const {
        return m_FactorCounts.size();
    }

    size_t CountOf(
        const mpz_class& Factor
    ) const {
        auto it = m_FactorCounts.find(Factor);
        if (it != m_FactorCounts.end()) {
            return it->second;
        }
        return 0;
    }

    bool Empty(
        void
    ) const {
        return m_FactorCounts.empty();
    }

    void Clear(
        void
    ) {
        m_FactorCounts.clear();
    }

    std::vector<std::pair<mpz_class, size_t>> ToVector() const {
        std::vector<std::pair<mpz_class, size_t>> vec;
        for (const auto& pair : m_FactorCounts) {
            vec.emplace_back(pair.first, pair.second);
        }
        return vec;
    }

    // Function to convert prime factors to a vector of composite factors
    // By iterating over all combinations of prime factor powers
    std::vector<mpz_class>
    GetComposite(
        const bool Sorted = false
    ) const {
        std::vector<mpz_class> composites;
        composites.push_back(1); // Start with 1 as the first composite
        for (const auto& [prime, count] : m_FactorCounts) {
            size_t current_size = composites.size();
            mpz_class prime_power = 1;
            for (size_t i = 1; i <= count; ++i) {
                prime_power *= prime;
                for (size_t j = 0; j < current_size; ++j) {
                    composites.push_back(composites[j] * prime_power);
                }
            }
        }
        if (Sorted)
        {
            std::sort(composites.begin(), composites.end());
        }
        return composites;
    }

    mpz_class
    Product(
        void
    ) const {
        mpz_class prod = 1;
        mpz_class factor;
        for (const auto& [prime, count] : m_FactorCounts) {
            mpz_pow_ui(factor.get_mpz_t(), prime.get_mpz_t(), count);
            prod *= factor;
        }
        return prod;
    }

    uint64_t
    Product64(
        void
    ) const {
        uint64_t prod = 1;
        mpz_class factor;
        for (const auto& [prime, count] : m_FactorCounts) {
            mpz_pow_ui(factor.get_mpz_t(), prime.get_mpz_t(), count);
            prod *= factor.get_ui();
        }
        return prod;
    }

    std::vector<uint8_t>
    Serialize(
        void
    ) const {
        std::vector<uint8_t> data;
        return data;
    }

    mpz_class
    MaxFactor(
        void
    ) const {
        if (m_FactorCounts.empty()) {
            return 0;
        }
        return std::max_element(
            m_FactorCounts.begin(),
            m_FactorCounts.end(),
            [](const auto& a, const auto& b) {
                return a.first < b.first;
            })->first;
    }

    size_t
    MaxPower(
        void
    ) const {
        if (m_FactorCounts.empty()) {
            return 0;
        }
        return std::max_element(
            m_FactorCounts.begin(),
            m_FactorCounts.end(),
            [](const auto& a, const auto& b) {
                return a.second < b.second;
            })->second;
    }
private:
    std::map<mpz_class, size_t> m_FactorCounts;
};