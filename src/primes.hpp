#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

#include <gmpxx.h>

const bool
load_prime_gaps(
    std::string_view filename
);

std::vector<uint8_t>
generate_prime_gaps(
    const mpz_class& limit,
    const bool is_count
);

std::optional<std::span<const uint8_t>>
get_prime_gaps(
    const uint64_t fallback_limit = 65536ull
);

mpz_class
get_nth_prime(
    const size_t n
);

size_t
get_prime_index(
    const mpz_class& prime
);