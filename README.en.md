# Wild Goppa Codes over Prime Fields

Experimental C implementation of a Goppa-code construction over the prime field
\\(\\mathbb{F}_N\\), with syndrome computation and an attempted
Berlekamp-style locator-polynomial decoder.

## Parameters

The current defaults are defined in `global-p.h`:

- `N`: prime-field modulus and code length
- `K`: polynomial/constraint parameter
- `T`: configured error weight (`K / 2`)
- `E`: auxiliary parameter used by the implementation

These values are experimental. In particular, changing `N` requires it to be
prime for the modular inverses used throughout the code to be well-defined.

## Build and run

On Linux with GCC and the standard C development toolchain:

```sh
gcc -O3 -Wall -Wextra -o decoder fujiyama.c -lm
ulimit -s unlimited
./decoder
```

The program uses large statically sized arrays and x86 intrinsics are included
in the source. A compatible x86 GCC environment may be required.

## Current implementation notes

- `fujiyama.c` contains the polynomial and finite-field arithmetic, parity-check
  matrix construction, syndrome calculation, locator-polynomial routines, and
  the demonstration `main`.
- `global-p.h` holds the compile-time parameters and shared matrices.
- `struct.h` defines the polynomial and matrix data structures.
- `chash-p.c` is included by `fujiyama.c` and contains supporting routines.

The demonstration now generates a weight-`T` error vector, runs the decoder,
and compares the recovered vector with the original. A successful comparison
prints `PASS`; a mismatch or failed Forney evaluation returns a failure status.

**Important:** this is research/experimental code, not a validated cryptographic
implementation. The decoding path and its algebraic conventions still require
independent mathematical review and reproducible testing across parameter sets.
Do not use it to protect data.

## Why prime-field Goppa codes?

- The prime modulus gives more flexibility in choosing the code length than
  restricting the length to powers of two.
- Arithmetic is expressed with integer operations modulo a prime, which can
  make small prototypes easier to inspect and experiment with.
