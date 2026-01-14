#include <string_view>
#include <iostream>

#include "primefactorcache.hpp"

int main(
    int argc,
    char* argv[]
)
{
    if (argc < 2) {
        std::cerr << "Usage: cachesort <cache_path>" << std::endl;
        return 1;
    }

    std::string_view cache_path = argv[1];
    
    // Create cache in read mode just to access the sort function
    PrimeFactorCache cache(cache_path);
    
    std::cout << "Sorting cache at: " << cache_path << std::endl;
    cache.Sort();
    std::cout << "Done sorting cache." << std::endl;
    
    return 0;
}
