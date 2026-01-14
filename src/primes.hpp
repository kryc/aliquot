#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

#include <gmpxx.h>

constexpr size_t kBitsPerWheelGap = 6;
constexpr uint64_t kGapMask = (1 << kBitsPerWheelGap) - 1;
constexpr size_t kMaxWheelGap = (1 << kBitsPerWheelGap) - 1;
constexpr size_t kGapsPerWord = (sizeof(uint64_t) * 8) / kBitsPerWheelGap;

constexpr uint32_t kWheel30 = 0x26424246;
constexpr uint32_t kWheel30BitsPerGap = 4;
constexpr uint32_t kWheel30Mask = (1 << kWheel30BitsPerGap) - 1;

const bool
load_prime_gaps(
    std::string_view filename
);

std::vector<uint8_t>
generate_prime_gaps(
    const mpz_class& limit,
    const bool is_count
);

std::span<const uint8_t>
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

std::span<const uint64_t>
get_primes_for_wheel_modulus(
    const size_t modulus
);

std::span<uint64_t>
get_wheel(
    const size_t modulus
);