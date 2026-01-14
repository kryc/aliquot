#include <array>
#include <cinttypes>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

#include <gmpxx.h>
#include <sys/mman.h>

#include "primes.hpp"

static std::vector<uint8_t> gGeneratedPrimeGaps;
static std::span<const uint8_t> gMappedPrimes;
static std::string_view gPrimesFilename;
static FILE* gPrimesFile = nullptr;

static const std::array<uint64_t, 12> gFirstPrimes = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};

static std::map<size_t, std::vector<uint64_t>> gWheelCache;

const bool
LoadPrimeGaps(
    std::string_view Filename
)
{
    if (gMappedPrimes.data() != nullptr) {
        // Unmap the previous file
        munmap(
            const_cast<uint8_t*>(gMappedPrimes.data()),
            gMappedPrimes.size()
        );
        gMappedPrimes = std::span<const uint8_t>();
        gPrimesFilename = "";
        std::fclose(gPrimesFile);
        gPrimesFile = nullptr;
    }

    // Mmap the file
    const size_t size = std::filesystem::file_size(Filename);
    FILE* file = std::fopen(Filename.data(), "rb");
    if (!file) {
        return false;
    }

    const uint8_t* base = static_cast<const uint8_t*>(
        mmap(
            nullptr,
            size,
            PROT_READ,
            MAP_PRIVATE,
            fileno(file),
            0
        )
    );
    if (base == MAP_FAILED) {
        fclose(file);
        return false;
    }
    gMappedPrimes = std::span<const uint8_t>(base, size);
    gPrimesFilename = Filename;
    gPrimesFile = file;
    return true;
}

std::vector<uint8_t>
GeneratePrimeGaps(
    const mpz_class& Limit,
    const bool IsCount
)
{
    std::vector<uint8_t> gaps;
    if (Limit < 3) {
        return gaps;
    }

    gaps.reserve(10 * 1024 * 1024); // Reserve 100MB initially

    mpz_class next = 2;
    mpz_class previous = 0;
    size_t count = 0;
    while ((count && count < Limit) ||
        (!count && next <= Limit)) {
        mpz_class gap = next - previous;
        previous = next;
        // VLE encode the gap
        while (gap > 0) {
            uint8_t byte = gap.get_ui() & 0x7F;
            gap >>= 7;
            if (gap > 0) {
                byte |= 0x80; // More bytes to come
            }
            gaps.push_back(byte);
        }
        count++;
        mpz_nextprime(next.get_mpz_t(), next.get_mpz_t());
    }
    return gaps;
}

std::span<const uint8_t>
GetPrimeGaps(
    const uint64_t FallbackLimit
)
{
    if (gMappedPrimes.empty() && gPrimesFilename.empty()) {
        gGeneratedPrimeGaps = GeneratePrimeGaps(FallbackLimit, false);
        return gGeneratedPrimeGaps;
    } else if (gMappedPrimes.empty()) {
        // Failed to mmap the file
        return gGeneratedPrimeGaps;
    }
    return gMappedPrimes;
}

mpz_class
GetNthPrime(
    const size_t N
)
{
    if (N == 0) {
        return 2;
    }
    
    const auto& gaps = GetPrimeGaps();

    mpz_class prime = 2;
    size_t gap_index = 1; // Start after the first gap which is for prime 2
    size_t count = 0;

    while (count < N && gap_index < gaps.size()) {
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
        count++;
    }

    // Continue with fallback to mpz_nextprime if needed
    // if (count < N) {
    //     std::cout << "Continuing to find prime " << N << " using mpz_nextprime. Currently at " << prime << std::endl;
    // }
    while (count < N) {
        mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
        count++;
    }

    return prime;
}

size_t
GetPrimeIndex(
    const mpz_class& Prime
)
{
    const auto gaps = GetPrimeGaps();

    mpz_class current_prime = 2;
    size_t gap_index = 1; // Start after the first gap which is for prime 2
    size_t index = 0;

    while (current_prime < Prime && gap_index < gaps.size()) {
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
        current_prime += gap;
        index++;
    }

    // Fallback using mpz_nextprime
    while (current_prime < Prime) {
        mpz_nextprime(current_prime.get_mpz_t(), current_prime.get_mpz_t());
        index++;
    }

    return index;
}

std::span<const uint64_t>
GetPrimesForWheelModulus(
    const size_t Modulus
)
{
    std::span<const uint64_t> primes = gFirstPrimes;
    switch(Modulus)
    {
        case 30:
            return primes.subspan(0, 3); // 2 - 5
        case 210:
            return primes.subspan(0, 4); // 2 - 7
        case 2310:
            return primes.subspan(0, 5); // 2 - 11
        case 30030:
            return primes.subspan(0, 6); // 2 - 13
        case 510510:
            return primes.subspan(0, 7); // 2 - 17
        case 9699690:
            return primes.subspan(0, 8); // 2 - 19
        case 223092870:
            return primes.subspan(0, 9); // 2 - 23
        case 6469693230:
            return primes.subspan(0, 10); // 2 - 29
        default:
            throw std::invalid_argument("Unsupported wheel Modulus");
    }
}

std::span<const uint64_t>
GetWheel(
    const size_t Modulus
)
{
    // Check cache first
    if (gWheelCache.find(Modulus) != gWheelCache.end()) {
        return gWheelCache[Modulus];
    }

    // Generate the gaps for the wheel
    std::vector<uint64_t> gaps;
    uint64_t next = 0;
    size_t count = 0;
    uint64_t last_residue = 1;
    for (size_t residue = 3; residue < Modulus; residue+=2) {
        if (std::gcd(residue, Modulus) != 1) {
            continue;
        }
        uint64_t gap = residue - last_residue;
        if (gap > kMaxWheelGap) {
            throw std::runtime_error("Wheel gap exceeds maximum representable size.");
        }
        const size_t shift = count++ * kBitsPerWheelGap;
        next |= (gap << shift);
        if (count == kGapsPerWord) {
            gaps.push_back(next);
            next = 0;
            count = 0;
        }
        last_residue = residue;
    }

    // Add final gap to complete the cycle
    const uint64_t gap = Modulus - last_residue + 1;
    if (count > 0) {
        const size_t shift = count++ * kBitsPerWheelGap;
        next |= (gap << shift);
    } else {
        next = gap;
    }

    gaps.push_back(next);
    // Cache the result and return
    gWheelCache[Modulus] = gaps;
    return gWheelCache[Modulus];
}