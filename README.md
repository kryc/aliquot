# Aliqout

A C++ application for generating [aliquot sequences](https://en.wikipedia.org/wiki/Aliquot_sequence).

# Building

## Requirements

```bash
sudo apt install build-essential clang-21 cmake libssl-dev libicu-dev libgmp-dev
```

## Building

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

# Usage

## aliquot
```bash
./aliquot --help
Usage: aliquot [options] <number>
Options:
  -p <file>   Load prime gaps from file
  -h, --help  Show this help message
```

For example, the sequece for 24:

```bash
./aliquot 24
Aliquot sequence for 24:
36
55
17
1
```

## primegen

For larger numbers, we can use pre-computed prime numbers to speed up the prime factorisation step.

The project contains a binary called `primegen` for this purpose.

```bash
./primegen --help
Usage: primegen [options] <output_file>
Options:
  -2 <N>    Generate primes up to 2^N
  -n <N>    Generate primes up to N
  -c <N>    Generate first N primes
```

For example, to generate all primes up to 2<sup>36</sup>

```bash
./primegen -2 36
```

They are stored extremely efficiently using variable-length encoding, with an average of one byte per prime number. The above file uses less than 3GB of disk space.

Using `primegen` is optional but greatly speeds up processing.