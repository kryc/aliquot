
#include <filesystem>
#include <fstream>
#include <iostream>
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

const bool
load_prime_gaps(
    std::string_view filename
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
    const size_t size = std::filesystem::file_size(filename);
    FILE* file = std::fopen(filename.data(), "rb");
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
    gPrimesFilename = filename;
    gPrimesFile = file;
    return true;
}

std::vector<uint8_t>
generate_prime_gaps(
    const mpz_class& limit,
    const bool is_count
)
{
    std::vector<uint8_t> gaps;
    if (limit < 3) {
        return gaps;
    }

    gaps.reserve(10 * 1024 * 1024); // Reserve 100MB initially

    mpz_class next = 2;
    mpz_class previous = 0;
    size_t count = 0;
    while ((count && count < limit) ||
        (!count && next <= limit)) {
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
get_prime_gaps(
    const uint64_t fallback_limit
)
{
    if (gMappedPrimes.empty() && gPrimesFilename.empty()) {
        gGeneratedPrimeGaps = generate_prime_gaps(fallback_limit, false);
        return gGeneratedPrimeGaps;
    } else if (gMappedPrimes.empty()) {
        // Failed to mmap the file
        return gGeneratedPrimeGaps;
    }
    return gMappedPrimes;
}

mpz_class
get_nth_prime(
    const size_t n
)
{
    if (n == 0) {
        return 2;
    }
    
    const auto& gaps = get_prime_gaps();

    mpz_class prime = 2;
    size_t gap_index = 1; // Start after the first gap which is for prime 2
    size_t count = 0;

    while (count < n && gap_index < gaps.size()) {
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
    // if (count < n) {
    //     std::cout << "Continuing to find prime " << n << " using mpz_nextprime. Currently at " << prime << std::endl;
    // }
    while (count < n) {
        mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
        count++;
    }

    return prime;
}

size_t
get_prime_index(
    const mpz_class& prime
)
{
    const auto gaps = get_prime_gaps();

    mpz_class current_prime = 2;
    size_t gap_index = 1; // Start after the first gap which is for prime 2
    size_t index = 0;

    while (current_prime < prime && gap_index < gaps.size()) {
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
    while (current_prime < prime) {
        mpz_nextprime(current_prime.get_mpz_t(), current_prime.get_mpz_t());
        index++;
    }

    return index;
}