#pragma once

#include <algorithm>
#include <cinttypes>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <span>
#include <sys/mman.h>
#include <string_view>
#include <vector>

#include "factors.hpp"

typedef struct _index_entry {
    __uint128_t product;
    size_t num_factors;
} IndexEntry;

typedef struct _factor {
    __uint128_t value;
    size_t count;
} Factor;

struct FactorRecord {
    __uint128_t product;
    Factor factors[0];  // Flexible array member
};

class PrimeFactorCache {
public:
    PrimeFactorCache(
        const std::string_view path = ""
    ) : cache_path(path) {
        if (cache_path.empty()) {
            return;
        }

        if (!std::filesystem::exists(cache_path)) {
            std::filesystem::create_directories(cache_path);
        }

        if (!std::filesystem::exists(get_index_path())) {
            std::filesystem::create_directories(get_index_path());
        }
    };

    ~PrimeFactorCache() {
        close();
    }

    const bool
    is_open(
        void
    ) const {
        return !cache_path.empty();
    }

    std::filesystem::path
    get_path(
        void
    ) const {
        return cache_path;
    }

    std::filesystem::path
    get_index_path(
        void
    ) const {
        return cache_path / "index";
    }

    std::filesystem::path
    get_index_path(
        const uint8_t lowbyte
    ) const {
        return cache_path / "index" / (std::to_string(lowbyte) + ".idx");
    }

    std::filesystem::path
    get_factor_path(
        const size_t num_factors
    ) const {
        return cache_path / ("factors_" + std::to_string(num_factors) + ".dat");
    }

    std::filesystem::path
    get_info_path(
        void
    ) const {
        return cache_path / "info.txt";
    }

    void
    write_info(
        const size_t min_prime,
        const size_t max_prime,
        const size_t min_factors,
        const size_t max_factors,
        const size_t smallest_factor
    ) const {
        std::ofstream info_file(get_info_path(), std::ios::trunc);
        if (info_file.is_open()) {
            info_file << "min_prime=" << min_prime << "\n";
            info_file << "max_prime=" << max_prime << "\n";
            info_file << "min_factors=" << min_factors << "\n";
            info_file << "max_factors=" << max_factors << "\n";
            info_file << "smallest_factor=" << smallest_factor << "\n";
            info_file.close();
        }
    }

    std::optional<PrimeFactors>
    product_exists(
        const uint64_t product
    ) {
        // Get the low bytes to find the correct index span
        uint8_t lowbyte = static_cast<uint8_t>(product & 0xFF);
        // Check if the index exists
        if (std::filesystem::exists(get_index_path(lowbyte)) == false) {
            return std::nullopt;
        }
        // Get the file size
        const size_t index_file_size = std::filesystem::file_size(get_index_path(lowbyte));
        const size_t num_entries = index_file_size / sizeof(IndexEntry);
        if (num_entries == 0) {
            return std::nullopt;
        }
        // Open a handle to the index file
        FILE* indexfd = fopen(get_index_path(lowbyte).c_str(), "r");
        if (indexfd == nullptr) {
            throw std::runtime_error("Failed to open index file for reading: " + get_index_path(lowbyte).string());
        }

        // Binary search in the index
        ssize_t low = 0;
        ssize_t high = static_cast<ssize_t>(num_entries) - 1;
        size_t num_factors = 0;
        while (low <= high) {
            ssize_t mid = low + (high - low) / 2;
            // Seek to the record position
            if (fseek(indexfd, mid * sizeof(IndexEntry), SEEK_SET) != 0) {
                throw std::runtime_error("Failed to seek in index file.");
            }
            IndexEntry entry;
            // Read the record
            if (fread(&entry, sizeof(IndexEntry), 1, indexfd) != 1) {
                throw std::runtime_error("Failed to read record from index file.");
            }
            if (entry.product == product) {
                num_factors = entry.num_factors;
                break;
            } else if (entry.product < product) {
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
        fclose(indexfd);
    
        if (num_factors == 0) {
            return std::nullopt;
        }

        // Read the factor record from the appropriate factor file
        std::filesystem::path factor_path = get_factor_path(num_factors);
        FILE* factor_fd = fopen(factor_path.c_str(), "r");
        if (factor_fd == nullptr) {
            throw std::runtime_error("Failed to open factor file for reading: " + factor_path.string());
        }

        // Get file size and calculate record parameters
        const size_t factor_file_size = std::filesystem::file_size(get_factor_path(num_factors));
        const size_t record_size = sizeof(FactorRecord) + num_factors * sizeof(Factor);
        size_t num_records = factor_file_size / record_size;
        
        // Allocate buffer for reading one record at a time
        std::vector<char> record_buffer(record_size);
        FactorRecord* record = reinterpret_cast<FactorRecord*>(record_buffer.data());
        
        ssize_t factor_low = 0;
        ssize_t factor_high = static_cast<ssize_t>(num_records) - 1;
        
        while (factor_low <= factor_high) {
            ssize_t mid = factor_low + (factor_high - factor_low) / 2;
            
            // Seek to the record position
            if (fseek(factor_fd, mid * record_size, SEEK_SET) != 0) {
                throw std::runtime_error("Failed to seek in factor file.");
            }
            
            // Read the record
            if (fread(record_buffer.data(), record_size, 1, factor_fd) != 1) {
                throw std::runtime_error("Failed to read record from factor file.");
            }
            
            if (record->product == product) {
                PrimeFactors factors;
                for (size_t i = 0; i < num_factors; ++i) {
                    mpz_class prime_value(record->factors[i].value);
                    for (size_t j = 0; j < record->factors[i].count; ++j) {
                        factors.add_factor(prime_value);
                    }
                }
                fclose(factor_fd);
                return factors;
            } else if (record->product < product) {
                factor_low = mid + 1;
            } else {
                factor_high = mid - 1;
            }
        }
        fclose(factor_fd);
        return std::nullopt;
    }

    void write(
        const PrimeFactors factors
    ) {
        const size_t num_factors = factors.size();
        const uint64_t product = factors.product64();
        const uint8_t lowbyte = static_cast<uint8_t>(product & 0xFF);
        // Open a handle to the index file
        FILE* indexfd = fopen(get_index_path(lowbyte).c_str(), "a");
        if (indexfd == nullptr) {
            throw std::runtime_error("Failed to open index file for writing: " + get_index_path(lowbyte).string());
        }
        // Write the index entry
        IndexEntry entry;
        entry.product = product;
        entry.num_factors = static_cast<uint8_t>(num_factors);
        if (fwrite(&entry, sizeof(IndexEntry), 1, indexfd) != 1) {
            throw std::runtime_error("Failed to write index entry.");
        }
        fclose(indexfd);
        // Sort this index file
        sort_index(lowbyte);
        
        // Allocate buffer for FactorRecord with flexible array member
        const size_t record_size = sizeof(FactorRecord) + num_factors * sizeof(Factor);
        std::vector<char> buffer(record_size);
        FactorRecord* record = reinterpret_cast<FactorRecord*>(buffer.data());
        
        record->product = product;
        size_t i = 0;
        for (const auto& [prime, count] : factors.to_vector()) {
            record->factors[i].value = prime.get_ui();
            record->factors[i].count = static_cast<uint8_t>(count);
            i++;
        }
        
        // Write to the appropriate factor file
        std::filesystem::path factor_path = get_factor_path(num_factors);
        FILE* factor_fd = fopen(factor_path.c_str(), "a");
        if (factor_fd == nullptr) {
            throw std::runtime_error("Failed to open factor file for writing: " + factor_path.string());
        }
        if (fwrite(buffer.data(), record_size, 1, factor_fd) != 1) {
            throw std::runtime_error("Failed to write factor record.");
        }
        fclose(factor_fd);
        // Sort this factor file
        sort_factors(num_factors);
    }

    void close(
        void
    ) {}

    bool sort_index(
        const size_t lowbyte
    ) const {
        std::filesystem::path index_path = get_index_path(static_cast<uint8_t>(lowbyte));
        if (!std::filesystem::exists(index_path)) {
            return false;
        }
        size_t file_size = std::filesystem::file_size(index_path);
        if (file_size == 0) {
            return false;  // Skip empty index files
        }
        // Open and mmap the index file
        FILE* fd = fopen(index_path.c_str(), "r+");
        if (fd == nullptr) {
            std::cerr << "Failed to open index file for sorting." << std::endl;
            return false;
        }
        void* map = mmap(nullptr, file_size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(fd), 0);
        if (map == MAP_FAILED) {
            std::cerr << "Failed to mmap index file for sorting." << std::endl;
            fclose(fd);
            return false;
        }
        size_t num_entries = file_size / sizeof(IndexEntry);
        std::span<IndexEntry> entries(static_cast<IndexEntry*>(map), num_entries);
        // Sort the entries by product
        std::sort(entries.begin(), entries.end(), [](const IndexEntry& a, const IndexEntry& b) {
            return a.product < b.product;
        });
        // Sync changes to disk
        msync(map, file_size, MS_SYNC);
        // Unmap and close the file
        munmap(map, file_size);
        fclose(fd);
        return true;
    }

    bool sort_factors(
        const size_t num_factors
    ) const {
        std::filesystem::path factor_path = get_factor_path(num_factors);
        if (!std::filesystem::exists(factor_path)) {
            return false;
        }
        // Open and mmap the factor file
        FILE* ffd = fopen(factor_path.c_str(), "r+");
        if (ffd == nullptr) {
            std::cerr << "Failed to open factor file for sorting: " << factor_path << std::endl;
            return false;
        }
        size_t factor_file_size = std::filesystem::file_size(factor_path);
        void* fmap = mmap(nullptr, factor_file_size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(ffd), 0);
        if (fmap == MAP_FAILED) {
            std::cerr << "Failed to mmap factor file for sorting: " << factor_path << std::endl;
            fclose(ffd);
            return false;
        }
        
        const size_t record_size = sizeof(FactorRecord) + num_factors * sizeof(Factor);
        size_t num_records = factor_file_size / record_size;
        
        // Comparison function for qsort_r
        auto compare_records = [](const void* a, const void* b, void* arg) -> int {
            const FactorRecord* rec_a = static_cast<const FactorRecord*>(a);
            const FactorRecord* rec_b = static_cast<const FactorRecord*>(b);
            if (rec_a->product < rec_b->product) return -1;
            if (rec_a->product > rec_b->product) return 1;
            return 0;
        };
        
        // Sort in-place using qsort_r
        qsort_r(fmap, num_records, record_size, compare_records, nullptr);
        
        // Sync changes to disk
        msync(fmap, factor_file_size, MS_SYNC);
        // Unmap and close the file
        munmap(fmap, factor_file_size);
        fclose(ffd);
        return true;
    }

    void sort(
        void
    ) {
        // Sort index files
        for (size_t i = 0; i < 256; ++i) {
            sort_index(i);
        }

        for (size_t num_factors = 1; num_factors <= 6; ++num_factors) {
            sort_factors(num_factors);
        }
    }

    void print_stats(
        void
    ) const {
        std::cout << "Prime Factor Cache Stats:" << std::endl;
        std::cout << "Cache Path: " << cache_path << std::endl;
        // Get number of entries in the index file
        size_t index_size = 0;
        for (size_t i = 0; i < 256; ++i) {
            if (std::filesystem::exists(get_index_path(i))) {
                index_size = std::filesystem::file_size(get_index_path(i));
                index_size = index_size / sizeof(IndexEntry);
            }
        }
        std::cout << "Entries: " << index_size << std::endl;
        for (size_t i = 1; i <= 6; ++i) {
            std::filesystem::path factor_path = get_factor_path(i);
            if (std::filesystem::exists(factor_path)) {
                size_t factor_size = std::filesystem::file_size(factor_path);
                size_t record_size = sizeof(FactorRecord) + i * sizeof(Factor);
                size_t num_records = factor_size / record_size;
                std::cout << "Factors with " << i << " primes: " << num_records << std::endl;
            }
        }
    }
private:
    std::filesystem::path cache_path;
};