#ifndef RUTILS_H
#define RUTILS_H

#include <stdio.h>
#include <stdint.h>
#include "ru_ksort.h"

#ifdef __cplusplus
extern "C" {
#endif

#define sort_key_128x(a) ((a).x)
#define sort_key_64(x) (x)

extern double ri_realtime0;
extern int ri_verbose;
extern unsigned char seq_nt4_table[256];

// emulate 128-bit integers and arrays
typedef struct { uint64_t x, y; } mm128_t;
typedef struct { uint32_t x, y; } mm64_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;
typedef struct { size_t n, m; mm64_t *a; } mm64_v;
typedef struct { size_t n, m; char **a; } ri_char_v;

void radix_sort_128x(mm128_t *beg, mm128_t *end);
void radix_sort_64(uint64_t *beg, uint64_t *end);
void radix_sort_64x(mm64_v *beg, mm64_v *end);
uint32_t ks_ksmall_uint32_t(size_t n, uint32_t arr[], size_t kk);

double ri_realtime(void);
double ri_cputime(void);
long ri_peakrss(void);

void load_pore(const char* fpore, const short k, const short lev_col, float** pore_vals);

#ifdef __cplusplus
}
#endif
#endif //RUTILS_H