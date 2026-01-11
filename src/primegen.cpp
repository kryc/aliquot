// primegen
// Precomputes the gaps between primes and outputs them to a file.
// The gaps are encoded as VLE-encoded bytes.
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string_view>
#include <cstdint>

#include <gmpxx.h>

std::string human_readable_size(size_t bytes) {
    const char* suffixes[] = {"B", "KB", "MB", "GB", "TB"};
    size_t s = 0;
    double count = static_cast<double>(bytes);
    while (count >= 1024 && s < 4) {
        s++;
        count /= 1024;
    }
    char buffer[50];
    snprintf(buffer, sizeof(buffer), "%.2f %s", count, suffixes[s]);
    return std::string(buffer);
}

int main(int argc, char* argv[]) {
    
    if (argc < 3) {
        std::cerr << "Usage: primegen [options] <output_file>" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << "  -2 <N>    Generate primes up to 2^N" << std::endl;
        std::cerr << "  -n <N>    Generate primes up to N" << std::endl;
        std::cerr << "  -c <N>    Generate first N primes" << std::endl;
        return 1;
    }

    std::string_view output_file;
    mpz_class max_prime = 0;
    bool use_count = false;

    for (int i = 1; i < argc; ++i) {
        if (std::string_view(argv[i]) == "-2" && i + 1 < argc) {
            mpz_class exp(argv[++i]);
            mpz_class two = 2;
            mpz_pow_ui(max_prime.get_mpz_t(), two.get_mpz_t(), exp.get_ui());
        } else if (std::string_view(argv[i]) == "-n" && i + 1 < argc) {
            max_prime = mpz_class(argv[++i]);
        } else if (std::string_view(argv[i]) == "-c" && i + 1 < argc) {
            max_prime = mpz_class(argv[++i]);
            use_count = true;
        } else {
            output_file = argv[i];
            break;
        }
    }

    mpz_class value = 2;
    mpz_class last = 0;
    size_t count = 1;
    size_t filesize = 0;

    std::ofstream ofs(output_file.data(), std::ios::binary);
    if (!ofs) {
        std::cerr << "Error: Could not open output file." << std::endl;
        return 1;
    }

    while (
        (use_count && count < max_prime) ||
        (!use_count && value <= max_prime)) {
        mpz_class gap = value - last;
        last = value;

        // VLE encode the gap
        while (gap > 0) {
            uint8_t byte = gap.get_ui() & 0x7F;
            gap >>= 7;
            if (gap > 0) {
                byte |= 0x80; // More bytes to come
            }
            ofs.put(static_cast<char>(byte));
            filesize++;
        }
        count++;
        
        if (count % 100000 == 0) {
            // Calculate percentage
            double percent = 0.0;
            if (!use_count) {
                percent = (value.get_d() / max_prime.get_d()) * 100.0;
            } else {
                percent = (static_cast<double>(count) / max_prime.get_d()) * 100.0;
            }
            std::cerr << "\r#: " << count
                    << " (" << human_readable_size(filesize) << "), "
                    << "latest prime: " << value
                    << std::fixed << std::setprecision(2)
                    << " (" << percent << "%)" << std::flush;
        }
        mpz_nextprime(value.get_mpz_t(), value.get_mpz_t());
    }

    std::cerr << std::endl << "Finished generating primes." << std::endl;
    std::cerr << "Output file size: " << human_readable_size(filesize) << std::endl;
    return 0;
}