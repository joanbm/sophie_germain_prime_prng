/* MIT License

Copyright (c) 2019 Joan Bruguera Micó

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
/// Generates an uniform sample, using a pseudorandom number generator
/// based on Sophie-Germain safe primes
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>

/**************************
 * CONFIGURATION / LIMITS *
 **************************/

//#define TEST_VALUES

#ifndef TEST_VALUES
/// Those values must be picked correctly, so no overflow happens during the algorithm
/// Some basic checks are asserted at the beginning of main(), which can be
/// used as a reference, should those values be changes
#define num_t uint64_t
#define NUM_MAX UINT64_MAX
#define PRInum PRIu64
#define SCNnum SCNu64
#define bignum_t unsigned __int128

#define NUM_OBSERVATIONS_MAX ((num_t)UINT32_MAX)
#define SEED_MAX ((num_t)UINT16_MAX)

#define NUM_DIGITS_PER_OBSERVATION 15

/// This is the maximum distance between two Sophie-Germain safe primes
/// (similar to the 'maximal prime gap' concept but with Sophie-Germain safe primes)
/// which are less than 2*(NUM_OBSERVATIONS_MAX + SEED_MAX).

/// This guarantees that there is at least one Sophie-Germain safe between
/// NUM_OBSERVATIONS_MAX + s     * NUM_PRIME_GERMAIN_GAP_MAX and
/// NUM_OBSERVATIONS_MAX + (s+1) * NUM_PRIME_GERMAIN_GAP_MAX,
/// for 0 <= s <= SEED_MAX, so that for any admissible value of s,
/// we can generate a single Sophie-Germain safe prime, and thus
/// generate a different pseudorandom sequence
#define NUM_PRIME_GERMAIN_GAP_MAX ((num_t)17904)
#else
#define num_t uint16_t
#define NUM_MAX UINT16_MAX
#define PRInum PRIu16
#define SCNnum SCNu16
#define bignum_t uint32_t

#define NUM_OBSERVATIONS_MAX ((num_t)255)
#define SEED_MAX ((num_t)15)

#define NUM_DIGITS_PER_OBSERVATION 2

#define NUM_PRIME_GERMAIN_GAP_MAX 616
#endif

/*********************
 * GENERIC UTILITIES *
 *********************/

/// Returns the number of elements in the specified array
#define ARRAY_SIZE(x) ((sizeof(x)/sizeof((x)[0])))

#define VALUE_STRINGIFY(a) MACRO_STRINGIFY(a)
#define MACRO_STRINGIFY(a) #a

/*********************
 * NUMERIC UTILITIES *
 *********************/

/// Returns -1, 0 or 1 according to whether a<b, a==b, or a>b (for bsearch)
static int compare_num(const void *a_ptr, const void *b_ptr) {
    num_t a = *((num_t *)a_ptr);
    num_t b = *((num_t *)b_ptr);

    if (a == b) {
        return 0;
    }
    return (a < b) ? -1 : 1;
}

/// Parses the specified string into a number.
static bool parse_num(const char *str, num_t *dest) {
    // Don't allow trailing noise. See https://stackoverflow.com/a/21888827
    char trailing_detect;
    return sscanf(str, "%" SCNnum "%c", dest, &trailing_detect) == 1;
}

/// Computes (x*y mod p) without multiplication potentially causing overflow
static num_t mul_mod(num_t x, num_t y, num_t p) {
    return (num_t)(((bignum_t)x * (bignum_t)y) % p);
}

/// Computes (x**y mod p) efficiently using modular exponentiation
/// See: https://en.wikipedia.org/wiki/Modular_exponentiation
static num_t pow_mod(num_t x, num_t y, num_t p) {
    if (p == 1) {
        return 0;
    }

    num_t result = 1;
    for (num_t curr_x = x % p, curr_y = y;
         curr_y > 0;
         curr_x = mul_mod(curr_x, curr_x, p), curr_y /= 2) {
        if (curr_y % 2 == 1) {
            result = mul_mod(result, curr_x, p);
        }
    }
    return result;
}

/*******************************************************************
 * IMPLEMENTATION OF THE RABIN-MILLER DETERMINISTIC PRIMALITY TEST *
 *******************************************************************/

/// List of Rabin-Miller witnesses that ensure (with 100% certainty) that
/// the Rabin-Miller test works for all 64-bit unsigned integer inputs
/// See https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Testing_against_small_sets_of_bases
static const num_t rm_witnesses[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
};

/// Tests the specified candidate passes the Rabin-Miller primality test for a witness, where:
/// p_candidate is an odd integer > 3 to be tested for primality
/// d and r are such that 2^r*d = p_candidate - 1, with d odd.
/// witness is the witness to be checked as a witness for the primality test
/// See: https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test
static bool test_rm_witness(num_t p_candidate, num_t d, num_t r, num_t witness) {
    num_t x = pow_mod(witness, d, p_candidate);
    if (x == 1 || x == p_candidate - 1) {
        return true;
    }

    for (num_t j = 0; j < (num_t)(r - 1); j++) {
        x = mul_mod(x, x, p_candidate);
        if (x == p_candidate - 1) {
            return true;
        }
    }

   return false;
}

/// Checks if a given number is a prime number (true) or not (false)
/// using the Rabin-Miller primality test
/// See: https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test
static bool rm_primality_test(num_t p_candidate) {
    if (p_candidate <= rm_witnesses[ARRAY_SIZE(rm_witnesses)-1]) {
        return bsearch(&p_candidate, rm_witnesses,
            ARRAY_SIZE(rm_witnesses), sizeof(num_t), compare_num) != NULL;
    }
    if (p_candidate % 2 == 0) {
        return false;
    }

    num_t d = (num_t)(p_candidate - 1);
    num_t r = 0;
    while (d % 2 == 0) {
        r++;
        d /= 2;
    }

    for (size_t i = 0; i < ARRAY_SIZE(rm_witnesses); i++) {
        if (!test_rm_witness(p_candidate, d, r, rm_witnesses[i])) {
            return false;
        }
    }

    return true;
}

/*********************************************************************************
 * IMPLEMENTATION OF THE SOPHIE-GERMAIN SAFE PRIME PSEUDORANDOM NUMBER GENERATOR *
 *********************************************************************************/

/// Checks if a given integer is a Sophie-Germain safe prime (aka. q)
/// Sophie-Germain prime condition (p where q=p*2+1)
/// See: https://en.wikipedia.org/wiki/Sophie_Germain_prime#Pseudorandom_number_generation
static bool is_sophie_germain_safe_prime(num_t q_candidate) {
    num_t p_candidate = (num_t)((q_candidate - 1) / 2);
    // Associated maximally periodic reciprocal condition for p
    num_t max_recip_test = p_candidate % 20;

    return (max_recip_test == 3 || max_recip_test == 9 || max_recip_test == 11) &&
           rm_primality_test(q_candidate) &&
           rm_primality_test(p_candidate);
}

/// Generate a Sophie-Germain safe prime greater or equal than the given
/// (inclusive) lower bound
/// See: https://en.wikipedia.org/wiki/Sophie_Germain_prime#Pseudorandom_number_generation
static num_t generate_sophie_germain_safe_prime(num_t lower_bound) {
    for (num_t q_candidate = lower_bound; q_candidate != NUM_MAX; q_candidate++) {
        if (is_sophie_germain_safe_prime(q_candidate)) {
            return q_candidate;
        }
    }

    return 0; // Not possible to find safe prime in acceptable range
}

/// Generates an uniform sample, using a pseudorandom number generator
/// based on Sophie-Germain safe primes
static void generate_uniform_sophie(num_t num_observations, num_t seed) {
    // Generate the required Sophie-Germain safe prime. If the program is correctly
    // configured, this generates a different value of q for every seed,
    // and it is greater than NUM_OBSERVATIONS_MAX * NUM_DIGITS_PER_OBSERVATION + 1,
    // so it will generate (at least) NUM_OBSERVATIONS_MAX * NUM_DIGITS_PER_OBSERVATION digits
    num_t min_q = (num_t)(NUM_OBSERVATIONS_MAX * NUM_DIGITS_PER_OBSERVATION + 1 + seed * NUM_PRIME_GERMAIN_GAP_MAX);
    fprintf(stderr, "Looking for a Sophie-Germain safe prime q >= %" PRInum "\n", min_q);

    num_t found_q = generate_sophie_germain_safe_prime(min_q);
    fprintf(stderr, "Found a Sophie-Germain safe prime q = %" PRInum "\n", found_q);
    assert(found_q != 0 && found_q <= min_q + NUM_PRIME_GERMAIN_GAP_MAX &&
        "Invalid configuration: Numeric overflow and/or incorrect NUM_PRIME_GERMAIN_GAP_MAX.");

    // Generate the decimal expansion of 1/q, that is, our random digits,
    // using a simple long-division based decimal digit extraction algorithm
    fprintf(stderr, "Generating the decimal expansion of 1/%" PRInum "...\n", found_q);

    num_t r = 1;
    char observation[NUM_DIGITS_PER_OBSERVATION+3];
    sprintf(observation, "0.%0" VALUE_STRINGIFY(NUM_DIGITS_PER_OBSERVATION) PRInum, (num_t)0);

    for (num_t i = 0; i < num_observations; i++) {
        for (size_t j = 0; j < NUM_DIGITS_PER_OBSERVATION; j++) {
            observation[j+2] = (char)('0' + ((r * 10) / found_q));
            r = (num_t)((r * 10) % found_q);
        }
        puts(observation);
    }
}

/// Entry point of the application. Parses the command line inputs and calls the generator
int main(int argc, char *argv[]) {
    // Do some basic sanity checks on the configuration
    assert(sizeof(num_t) * 2 <= sizeof(bignum_t) && "Invalid configuration: bignum_t misconfigured");

    assert((SEED_MAX * NUM_PRIME_GERMAIN_GAP_MAX) / NUM_PRIME_GERMAIN_GAP_MAX == SEED_MAX &&
        "Invalid configuration: (SEED_MAX * NUM_PRIME_GERMAIN_GAP_MAX) overflows.");

    assert((SEED_MAX * NUM_PRIME_GERMAIN_GAP_MAX) / NUM_PRIME_GERMAIN_GAP_MAX == SEED_MAX &&
        "Invalid configuration: (SEED_MAX * NUM_PRIME_GERMAIN_GAP_MAX) overflows.");

    assert((NUM_OBSERVATIONS_MAX * NUM_DIGITS_PER_OBSERVATION) / NUM_DIGITS_PER_OBSERVATION == NUM_OBSERVATIONS_MAX &&
        "Invalid configuration: (NUM_OBSERVATIONS_MAX * NUM_DIGITS_PER_OBSERVATION) overflows.");

    assert(SEED_MAX * NUM_PRIME_GERMAIN_GAP_MAX + NUM_OBSERVATIONS_MAX * NUM_DIGITS_PER_OBSERVATION+ 1 >
           SEED_MAX * NUM_PRIME_GERMAIN_GAP_MAX && "Invalid configuration: "
           "(SEED_MAX * NUM_PRIME_GERMAIN_GAP_MAX + NUM_OBSERVATIONS_MAX * NUM_DIGITS_PER_OBSERVATION + 1) overflows.");

    // Print title, check and validate command line arguments
    fprintf(stderr, "PRNG Based on Sophie-Germain primes\n");
    fprintf(stderr, "-----------------------------------\n");

    num_t num_observations, seed;
    if (argc != 3 ||
        !parse_num(argv[1], &num_observations) || num_observations > NUM_OBSERVATIONS_MAX ||
        !parse_num(argv[2], &seed) || seed > SEED_MAX)
    {
        fprintf(stderr, "Usage: %s num_observations seed\n", argv[0]);
        fprintf(stderr, "    (where num_observations <= %" PRInum ")\n", NUM_OBSERVATIONS_MAX);
        fprintf(stderr, "    (where seed <= %" PRInum ")\n", SEED_MAX);
        return EXIT_FAILURE;
    }

    // Once we have a valid parametrization, run the core algorithm
    generate_uniform_sophie(num_observations, seed);
    return EXIT_SUCCESS;
}
