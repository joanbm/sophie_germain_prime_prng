#!/usr/bin/env make

sophie_germain_prime_prng: sophie_germain_prime_prng.c
	gcc -Ofast -Wall -Wextra -Wconversion -std=c99 sophie_germain_prime_prng.c -osophie_germain_prime_prng

clean:
	rm -f sophie_germain_prime_prng