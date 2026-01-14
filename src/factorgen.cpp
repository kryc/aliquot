#include <iostream>
#include <cmath>
#include <string_view>

#include "factors.hpp"
#include "primefactorcache.hpp"
#include "primefactors.hpp"

const size_t
GenerateFactors(
    const uint64_t MinPrime,
    const uint64_t MaxPrime,
    const size_t MinNumFactors,
    const size_t MaxNumFactors,
    const uint64_t SmallestFactor,
    const std::string_view Output
)
{
    // First generate a list of primes in the range [min_prime, max_prime]
    std::vector<uint64_t> primes;
    mpz_class candidate = MinPrime;
    while (candidate <= MaxPrime) {
        if (mpz_probab_prime_p(candidate.get_mpz_t(), 25)) {
            primes.push_back(candidate.get_ui());
        }
        mpz_nextprime(candidate.get_mpz_t(), candidate.get_mpz_t());
    }
    const size_t num_primes = primes.size();
    
    // Now find the index of the SmallestFactor
    size_t start_index = 0;
    for (size_t i = 0; i < primes.size(); ++i) {
        if (primes[i] >= SmallestFactor) {
            start_index = i;
            break;
        }
    }

    // Create the prime factor cache
    PrimeFactorCache cache(Output);

    // Calculte the total number of products we will generate
    size_t total_products = 0;
    for (size_t n = MinNumFactors; n <= MaxNumFactors; ++n) {
        total_products += static_cast<size_t>(std::pow(num_primes - start_index, n));
    }
    std::cerr << "Generating approximately " << total_products << " products..." << std::endl;

    // Now generate products of factors
    size_t calculated = 0;
    size_t num_factors = 2;
    PrimeFactors factors;
    if (MinNumFactors >= 2) {
        for (size_t a = start_index; a < num_primes; ++a) {
            for (size_t b = 0; b < num_primes; ++b) {
                factors.Clear();
                factors.AddFactor(primes[a]);
                factors.AddFactor(primes[b]);
                cache.Write(factors);
                calculated++;

                if (calculated % 1000000 == 0) {
                    std::cerr << "\r" << num_factors << "/" << MaxNumFactors << " "
                    << calculated << "/" << total_products << "("
                    << (static_cast<double>(calculated) / total_products) * 100 << "%)" << std::flush;
                }
            }
        }
    }

    // Handle three factors
    if (MaxNumFactors >= 3) {
        num_factors = 3;
        for (size_t a = start_index; a < num_primes; ++a) {
            for (size_t b = 0; b < num_primes; ++b) {
                for (size_t c = 0; c < num_primes; ++c) {
                    factors.Clear();
                    factors.AddFactor(primes[a]);
                    factors.AddFactor(primes[b]);
                    factors.AddFactor(primes[c]);
                    cache.Write(factors);
                    calculated++;

                    if (calculated % 1000000 == 0) {
                        std::cerr << "\r" << num_factors << "/" << MaxNumFactors << " "
                        << calculated << "/" << total_products << "("
                        << (static_cast<double>(calculated) / total_products) * 100 << "%)" << std::flush;
                    }
                }
            }
        }
    }

    // Handle four factors
    if (MaxNumFactors >= 4) {
        num_factors = 4;
        for (size_t a = start_index; a < num_primes; ++a) {
            for (size_t b = 0; b < num_primes; ++b) {
                for (size_t c = 0; c < num_primes; ++c) {
                    for (size_t d = 0; d < num_primes; ++d) {
                        factors.Clear();
                        factors.AddFactor(primes[a]);
                        factors.AddFactor(primes[b]);
                        factors.AddFactor(primes[c]);
                        factors.AddFactor(primes[d]);
                        cache.Write(factors);
                        calculated++;

                        if (calculated % 1000000 == 0) {
                            std::cerr << "\r" << num_factors << "/" << MaxNumFactors << " "
                            << calculated << "/" << total_products << "("
                            << (static_cast<double>(calculated) / total_products) * 100 << "%)" << std::flush;
                        }
                    }
                }
            }
        }
    }
    
    cache.Close();
    std::cerr << std::endl;
    std::cerr << "Generated " << calculated << " products." << std::endl;
    std::cerr << "Sorting cache files..." << std::endl;
    cache.Sort();
    std::cerr << "Done sorting cache files." << std::endl;
    cache.PrintStats();

    // Write the info
    cache.WriteInfo(
        MinPrime,
        MaxPrime,
        MinNumFactors,
        MaxNumFactors,
        SmallestFactor
    );

    return 0;
}

int main(
    int argc,
    char* argv[]
) {
    if (argc < 2) {
        std::cerr << "Usage: factorgen [options] <output>" << std::endl;
        return 1;
    }

    uint64_t min_prime = 7;
    uint64_t max_prime = 65537;
    uint64_t smallest_factor = 1007;
    size_t min_num_factors = 2;
    size_t max_num_factors = 4;
    std::string_view output;

    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        if (arg == "-m" && i + 1 < argc) {
            min_prime = static_cast<uint64_t>(std::stoull(argv[++i]));
        } else if (arg == "-M" && i + 1 < argc) {
            max_prime = static_cast<uint64_t>(std::stoull(argv[++i]));
        } else if (arg == "-f" && i + 1 < argc) {
            min_num_factors = static_cast<size_t>(std::stoul(argv[++i]));
        } else if (arg == "-F" && i + 1 < argc) {
            max_num_factors = static_cast<size_t>(std::stoul(argv[++i]));
        } else if (arg == "-2" && i + 1 < argc) {
            uint64_t power = static_cast<uint64_t>(std::stoul(argv[++i]));
            smallest_factor = uint64_t(1) << power;
        } else if (arg == "-n" && i + 1 < argc) {
            smallest_factor = static_cast<uint64_t>(std::stoull(argv[++i]));
        } else if (arg == "-h" || arg == "--help") {
            std::cout << "Usage: factorgen [options] <output>" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  -m <N>    Minimum prime to use for factorization (default 7)" << std::endl;
            std::cout << "  -M <N>    Maximum prime to use for factorization (default 65537)" << std::endl;
            std::cout << "  -f <N>    Minimum number of factors (default 2)" << std::endl;
            std::cout << "  -F <N>    Maximum number of factors (default 4)" << std::endl;
            std::cout << "  -2 <N>    Set smallest factor value to 2^N" << std::endl;
            std::cout << "  -n <N>    Set smallest factor value to N" << std::endl;
            return 0;
        } else {
            output = arg;
            break;
        }
    }

    if (max_prime == 0) {
        std::cerr << "Error: Maximum prime not specified." << std::endl;
        return 1;
    }

    if (output.empty()) {
        std::cerr << "Error: Output file not specified." << std::endl;
        return 1;
    }

    // Make sure min_prime and max_prime are prime
    if (min_prime % 2 == 0) {
        min_prime += 1;
    }
    // while (mpz_probab_prime_p(min_prime.get_mpz_t(), 25) == 0) {
    //     min_prime += 2;
    // }

    if (max_prime % 2 == 0) {
        max_prime -= 1;
    }
    // while (mpz_probab_prime_p(max_prime.get_mpz_t(), 25) == 0) {
    //     max_prime -= 2;
    // }

    // Ensure smallest_factor is between min_prime and max_prime
    if (smallest_factor < min_prime || smallest_factor > max_prime) {
        std::cerr << "Error: Smallest factor must be between min_prime and max_prime." << std::endl;
        return 1;
    }

    GenerateFactors(
        min_prime,
        max_prime,
        min_num_factors,
        max_num_factors,
        smallest_factor,
        output
    );

}