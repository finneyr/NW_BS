
/*

vi nw_bs.c ; gcc -O2 -Wall -Wextra -o nw_bs nw_bs.c

  Needleman–Wunsch global alignment with bisulfite-aware matching.

  Directed BS matching (reference base vs read base):
    CT mode: ref 'C' matches read 'C' or 'T'
    GA mode: ref 'G' matches read 'G' or 'A'
    N in either string is treated as match-any (you can change that)

  Tie-break in traceback: DIAG > UP > LEFT.

  Compile:
    gcc -O2 -Wall -Wextra -o nw_bs nw_bs.c

  Usage (simple):
    ./nw_bs --mode CT  REFSEQ READSEQ
    ./nw_bs --mode EITHER -m 2 -x -1 -g -2 REFSEQ READSEQ
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>

typedef enum { TB_STOP=0, TB_DIAG=1, TB_UP=2, TB_LEFT=3 } tb_dir_t;
typedef enum { BS_CT=0, BS_GA=1 } bs_mode_t;

typedef struct {
    int score;
    bs_mode_t mode_used;
    char *aln_ref;
    char *aln_read;
    char *mid;   /* '|' for match, ' ' otherwise */
} nw_result_t;

static void die(const char *msg) {
    fprintf(stderr, "ERROR: %s\n", msg);
    exit(1);
}

static inline int max3(int a, int b, int c) {
    int m = (a > b) ? a : b;
    return (m > c) ? m : c;
}

#define IDX(i,j,cols) ((i)*(cols) + (j))

static int *alloc_int_mat(size_t rows, size_t cols) {
    if (rows == 0 || cols == 0) return NULL;
    if (rows > SIZE_MAX / cols) return NULL;
    size_t n = rows * cols;
    return (int*)malloc(n * sizeof(int));
}

static tb_dir_t *alloc_dir_mat(size_t rows, size_t cols) {
    if (rows == 0 || cols == 0) return NULL;
    if (rows > SIZE_MAX / cols) return NULL;
    size_t n = rows * cols;
    return (tb_dir_t*)malloc(n * sizeof(tb_dir_t));
}

/* Directed bisulfite “match?”: reference base r vs read base q */
static inline int bs_match(char ref_base, char read_base, bs_mode_t mode)
{
    unsigned char r = (unsigned char)toupper((unsigned char)ref_base);
    unsigned char q = (unsigned char)toupper((unsigned char)read_base);

    /* treat N as wildcard (edit if you want strict behavior) */
    if (r == 'N' || q == 'N') return 1;

    if (mode == BS_CT) {
        if (r == 'C') return (q == 'C' || q == 'T');
        return (q == r);
    } else { /* BS_GA */
        if (r == 'G') return (q == 'G' || q == 'A');
        return (q == r);
    }
}

static inline int sub_score(char ref_base, char read_base,
                            bs_mode_t mode, int match_score, int mismatch_score)
{
    return bs_match(ref_base, read_base, mode) ? match_score : mismatch_score;
}

static nw_result_t needleman_wunsch_bs(const char *ref, const char *read,
                                       bs_mode_t mode,
                                       int match_score, int mismatch_score, int gap_score)
{
    if (!ref || !read) die("NULL sequence pointer");

    size_t n = strlen(ref);
    size_t m = strlen(read);
    size_t rows = n + 1;
    size_t cols = m + 1;

    int *dp = alloc_int_mat(rows, cols);
    tb_dir_t *tb = alloc_dir_mat(rows, cols);
    if (!dp || !tb) die("allocation failed");

    dp[IDX(0,0,cols)] = 0;
    tb[IDX(0,0,cols)] = TB_STOP;

    for (size_t i = 1; i <= n; i++) {
        dp[IDX(i,0,cols)] = (int)i * gap_score;
        tb[IDX(i,0,cols)] = TB_UP;
    }
    for (size_t j = 1; j <= m; j++) {
        dp[IDX(0,j,cols)] = (int)j * gap_score;
        tb[IDX(0,j,cols)] = TB_LEFT;
    }

    for (size_t i = 1; i <= n; i++) {
        for (size_t j = 1; j <= m; j++) {
            int ssub = sub_score(ref[i-1], read[j-1], mode, match_score, mismatch_score);

            int score_diag = dp[IDX(i-1,j-1,cols)] + ssub;
            int score_up   = dp[IDX(i-1,j,cols)]   + gap_score;
            int score_left = dp[IDX(i,j-1,cols)]   + gap_score;

            int best = max3(score_diag, score_up, score_left);
            dp[IDX(i,j,cols)] = best;

            /* tie-break: DIAG > UP > LEFT */
            tb_dir_t dir = TB_DIAG;
            if (best == score_diag) dir = TB_DIAG;
            else if (best == score_up) dir = TB_UP;
            else dir = TB_LEFT;

            tb[IDX(i,j,cols)] = dir;
        }
    }

    /* traceback */
    size_t cap = n + m + 1;
    char *a_ref  = (char*)malloc(cap);
    char *a_read = (char*)malloc(cap);
    char *mid    = (char*)malloc(cap);
    if (!a_ref || !a_read || !mid) die("allocation failed (traceback buffers)");

    size_t i = n, j = m, k = 0;
    while (!(i == 0 && j == 0)) {
        tb_dir_t dir = tb[IDX(i,j,cols)];
        if (dir == TB_DIAG) {
            char r = ref[i-1];
            char q = read[j-1];
            a_ref[k]  = r;
            a_read[k] = q;
            mid[k]    = bs_match(r, q, mode) ? '|' : ' ';
            i--; j--;
        } else if (dir == TB_UP) {
            a_ref[k]  = ref[i-1];
            a_read[k] = '-';
            mid[k]    = ' ';
            i--;
        } else { /* TB_LEFT */
            a_ref[k]  = '-';
            a_read[k] = read[j-1];
            mid[k]    = ' ';
            j--;
        }
        k++;
        if (k + 1 >= cap) die("traceback overflow (unexpected)");
    }

    a_ref[k] = a_read[k] = mid[k] = '\0';
    for (size_t p = 0; p < k/2; p++) {
        char t;
        t = a_ref[p];  a_ref[p]  = a_ref[k-1-p];  a_ref[k-1-p]  = t;
        t = a_read[p]; a_read[p] = a_read[k-1-p]; a_read[k-1-p] = t;
        t = mid[p];    mid[p]    = mid[k-1-p];    mid[k-1-p]    = t;
    }

    nw_result_t res;
    res.score = dp[IDX(n,m,cols)];
    res.mode_used = mode;
    res.aln_ref = a_ref;
    res.aln_read = a_read;
    res.mid = mid;

    free(dp);
    free(tb);
    return res;
}

static void free_result(nw_result_t *r) {
    if (!r) return;
    free(r->aln_ref);
    free(r->aln_read);
    free(r->mid);
    r->aln_ref = r->aln_read = r->mid = NULL;
}

static const char *mode_name(bs_mode_t m) {
    return (m == BS_CT) ? "CT" : "GA";
}

static void usage(const char *argv0) {
    fprintf(stderr,
        "Usage:\n"
        "  %s --mode CT|GA|EITHER [-m match] [-x mismatch] [-g gap] REFSEQ READSEQ\n\n"
        "Defaults: match=1 mismatch=-1 gap=-1\n"
        "Bisulfite directed rules:\n"
        "  CT: ref C matches read C/T\n"
        "  GA: ref G matches read G/A\n",
        argv0);
    exit(1);
}

int main(int argc, char **argv)
{
    int match = 1, mismatch = -1, gap = -1;
    int have_mode = 0;
    int either = 0;
    bs_mode_t mode = BS_CT;

    int i = 1;
    while (i < argc && argv[i][0] == '-') {
        if (strcmp(argv[i], "--mode") == 0 && i+1 < argc) {
            const char *m = argv[++i];
            if (strcmp(m, "CT") == 0) { mode = BS_CT; either = 0; }
            else if (strcmp(m, "GA") == 0) { mode = BS_GA; either = 0; }
            else if (strcmp(m, "EITHER") == 0) { either = 1; }
            else usage(argv[0]);
            have_mode = 1;
        } else if (strcmp(argv[i], "-m") == 0 && i+1 < argc) {
            match = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-x") == 0 && i+1 < argc) {
            mismatch = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-g") == 0 && i+1 < argc) {
            gap = atoi(argv[++i]);
        } else {
            usage(argv[0]);
        }
        i++;
    }

    if (!have_mode) usage(argv[0]);
    if (argc - i != 2) usage(argv[0]);

    const char *ref = argv[i];
    const char *read = argv[i+1];

    nw_result_t best;
    memset(&best, 0, sizeof(best));
    best.score = INT32_MIN;

    if (!either) {
        best = needleman_wunsch_bs(ref, read, mode, match, mismatch, gap);
    } else {
        nw_result_t r_ct = needleman_wunsch_bs(ref, read, BS_CT, match, mismatch, gap);
        nw_result_t r_ga = needleman_wunsch_bs(ref, read, BS_GA, match, mismatch, gap);

        /* pick higher score; tie-break prefer CT */
        if (r_ct.score >= r_ga.score) {
            best = r_ct;
            free_result(&r_ga);
        } else {
            best = r_ga;
            free_result(&r_ct);
        }
    }

    printf("mode:  %s\n", mode_name(best.mode_used));
    printf("score: %d\n", best.score);
    printf("%s\n", best.aln_ref);
    printf("%s\n", best.mid);
    printf("%s\n", best.aln_read);

    free_result(&best);
    return 0;
}

