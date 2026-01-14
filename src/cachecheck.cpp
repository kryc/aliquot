
#include <string_view>

#include "primefactorcache.hpp"

int main(
    int argc,
    char* argv[]
)
{
    if (argc < 3) {
        std::cerr << "Usage: cachecheck <cache_path> value" << std::endl;
        return 1;
    }

    std::string_view cache_path = argv[1];
    uint64_t value = static_cast<uint64_t>(std::stoull(argv[2]));
    
    PrimeFactorCache cache(cache_path);
    auto result = cache.ProductExists(value);
    if (result.has_value()) {
        std::cout << "Product " << value << " exists in cache." << std::endl;
        
        // Output the factors
        auto factors = result.value();
        std::cout << "Factors: ";
        bool first = true;
        for (const auto& [prime, count] : factors.ToVector()) {
            if (!first) std::cout << " * ";
            if (count == 1) {
                std::cout << prime;
            } else {
                std::cout << prime << "^" << count;
            }
            first = false;
        }
        std::cout << std::endl;
        
        // Verify the product
        uint64_t computed_product = factors.Product64();
        if (computed_product == value) {
            std::cout << "Verification: PASSED (product = " << computed_product << ")" << std::endl;
        } else {
            std::cout << "Verification: FAILED (expected " << value << ", got " << computed_product << ")" << std::endl;
            return 1;
        }
    } else {
        std::cout << "Product " << value << " does not exist in cache." << std::endl;
    }

    return 0;
}