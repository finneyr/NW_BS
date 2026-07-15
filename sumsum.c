

/*
vi sumsum.c ; gcc -Werror -Wall -O2 -o sumsum sumsum.c -lm

Reads numbers (one per line) from stdin.
Single-pass stats (no sorting):
  cnt, sum, mean, min, max, variance (sample + population), stdev, mean(|x|)
Uses:
  - uint64_t for big counts
  - long double for higher precision
  - strtold() for robust parsing
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <errno.h>
#include <math.h>

#define MAXLINE 4000
static char s[MAXLINE];

int main(void)
{
    uint64_t cnt = 0;

    long double sum = 0.0L;
    long double mean = 0.0L;   // running mean (Welford)
    long double M2 = 0.0L;     // sum of squares of deviations from mean
    long double mean_abs = 0.0L;
    uint64_t n_nan=0, n_inf=0;
    long double minv = 0.0L, maxv = 0.0L;
    int have_minmax = 0;

    while (fgets(s, MAXLINE, stdin)) {
        char *end = NULL;
        errno = 0;

        long double x = strtold(s, &end);
        if (!isfinite(x)) 
        {  // in C99+, isfinite is a macro and works with long double on most libs
              fprintf(stderr, "Warning: non-finite value ignored: %s", s);
              continue;
        }
        // Skip empty/whitespace/non-numeric lines
        if (end == s) continue;
        if (isnan(x)) { n_nan++; continue; }
        if (isinf(x)) { n_inf++; continue; }

        // Ignore out-of-range conversions
        if (errno == ERANGE) {
            fprintf(stderr, "Warning: out-of-range value ignored: %s", s);
            continue;
        }

        // Update count early (Welford uses new n)
        cnt++;

        // sum
        sum += x;

        // min/max
        if (!have_minmax) {
            minv = maxv = x;
            have_minmax = 1;
        } else {
            if (x < minv) minv = x;
            if (x > maxv) maxv = x;
        }

        // mean(|x|) streaming update
        {
            long double ax = fabsl(x);
            mean_abs += (ax - mean_abs) / (long double)cnt;
        }

        // Welford: update mean and M2 (single pass, numerically stable)
        {
            long double delta = x - mean;
            mean += delta / (long double)cnt;
            long double delta2 = x - mean;
            M2 += delta * delta2;
        }
    }

    if (cnt == 0) {
        fprintf(stderr, "No numeric values read.\n");
        return 1;
    }

    long double var_pop = M2 / (long double)cnt;
    long double var_samp = (cnt > 1) ? (M2 / (long double)(cnt - 1)) : 0.0L;

    long double sd_pop = sqrtl(var_pop);
    long double sd_samp = sqrtl(var_samp);

    // Use %.21Lg for a compact, high-precision print of long doubles
    printf("cnt:%" PRIu64 " ", cnt);
    printf("sum:%.21Lg ", sum);
    printf("mean:%.21Lg ", mean);
    printf("mean_abs:%.21Lg ", mean_abs);
    printf("min:%.21Lg ", minv);
    printf("max:%.21Lg ", maxv);
    printf("var_pop:%.21Lg  sd_pop:%.21Lg ", var_pop, sd_pop);
    printf("var_samp:%.21Lg sd_samp:%.21Lg\n", var_samp, sd_samp);

    return 0;
}

