# Pseudorandom number generator based on Sophie-Germain safe primes

A pseudorandom number generator based on Sophie-Germain safe primes. More information on how to generate pseudorandom integers through Sophie-Germain safe primes can be found on [the Wikipedia page on Sophie-Germain primes](https://en.wikipedia.org/wiki/Sophie_Germain_prime#Pseudorandom_number_generation).

This repository is an implementation of such a pseudorandom number generator. It generates a sample of a uniform distribution from 0 to 1 (i.e. of `U(0,1)`), though this can be easily adapted to, for instance, generate a random digit stream.

Though the core idea of the generator is simple, some details such as seeding (so that a seed can be specified to generate different pseudo-random sequences on different runs), or implementing the generator with respectable efficiency, are not trivial. The focus of the implementation is on achieving those, while behaving correctly for all accepted inputs (to my best effort). 

It is written in C99 (standard, except GCC's 128-bit integer extension is used) and licensed under the MIT license. The code is pretty well commented, so you are encouraged to read it.

Note that no warranty is implied on the quality of the generated output, and in particular, due to the fact that this is a pseudorandom number generator, **it should not be used for any cryptography purpose whatsoever**.

## Quickstart

GCC must be installed. To compile the program, type `make`. The usage can be seen by executing the generated executable (`sophie`). For example:

```
$ make
gcc -Ofast -Wall -Wextra -Wconversion -std=c99 sophie.c -osophie
$ ./sophie 20 12345
PRNG Based on Sophie-Germain primes
-----------------------------------
Looking for a Sophie-Germain safe prime q >= 4515992176
Found a Sophie-Germain safe prime q = 4515994607
Generating the decimal expansion of 1/4515994607...
0.583510471406548
0.338245170096590
0.418713934482783
0.141197847759077
0.272963625574143
0.649547542517691
0.939307490851483
0.206831030655391
0.568761454900470
0.832205314013121
0.893880950775014
0.953422418912113
0.638846070482032
0.353764500410086
0.517182592375048
0.865951845455780
0.013955184955782
0.211367020793255
0.787871670503082
0.157467244291596
```

Only the generated uniform sample is printed to `stdout`, so the output can be easily captured into a file or piped to another program.
