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
LoadPrimeGaps(
    std::string_view Filename
);

std::vector<uint8_t>
GeneratePrimeGaps(
    const mpz_class& Limit,
    const bool IsCount
);

std::span<const uint8_t>
GetPrimeGaps(
    const uint64_t FallbackLimit = 65536ull
);

mpz_class
GetNthPrime(
    const size_t N
);

size_t
GetPrimeIndex(
    const mpz_class& Prime
);

std::span<const uint64_t>
GetPrimesForWheelModulus(
    const size_t Modulus
);

std::span<const uint64_t>
GetWheel(
    const size_t Modulus
);