#include <iostream>

#include <gmpxx.h>

#include "aliquot.hpp"
#include "primes.hpp"

static const std::string_view HELP_STRING = R"(
Usage: aliquot [options] <number>
Options:
  -p <file>   Load prime gaps from file
  -h, --help  Show this help message
)";

int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << HELP_STRING << std::endl;
        return 1;
    }

    std::string_view prime_gaps;
    mpz_class number;

    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        if (arg == "-p" && i + 1 < argc) {
            prime_gaps = argv[++i];
            if (!load_prime_gaps(prime_gaps)) {
                std::cerr << "Failed to load prime gaps from " << prime_gaps << std::endl;
                return 1;
            }
        } else if (arg == "-h" || arg == "--help") {
            std::cout << HELP_STRING << std::endl;
            return 0;
        } else {
            number = argv[i];
        }
    }

    if (number == 0) {
        std::cerr << "Please provide a valid number greater than 0." << std::endl;
        return 1;
    }

    try {
        std::cout << "Aliquot sequence for " << number << ":" << std::endl;
        auto sequence = aliquot_sequence(number, true);
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Error during prime factorization: " << ex.what() << std::endl;
        return 1;
    }

}