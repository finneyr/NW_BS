# NW_BS

`nw_bs` is a small, self-contained C program for performing **bisulfite-aware global pairwise sequence alignment** using the Needleman–Wunsch algorithm.

The program is intended for a restricted use case: aligning one predefined reference sequence against one predefined query or read sequence while treating expected bisulfite-conversion differences as compatible matches.

It is **not** intended to replace genome-scale bisulfite aligners or optimized general-purpose sequence-alignment libraries.

## Scope

`nw_bs` is designed for situations in which:

* the reference and query sequences have already been selected;
* only a small number of pairwise alignments are required;
* deterministic global alignment is appropriate;
* bisulfite-aware matching rules must be explicit;
* the output should be easy to inspect and parse;
* a single portable C source file is preferable to a larger software dependency.

The program does not:

* map FASTQ reads to a genome;
* search for candidate genomic locations;
* build or use genome indexes;
* perform methylation calling;
* process BAM or SAM files;
* support local alignment;
* implement affine-gap scoring;
* replace Bismark, bwa-meth, BSMAP, or similar genome-scale tools.

## Bisulfite-aware matching

The matching rule is directed from the unconverted reference sequence to the potentially bisulfite-converted query sequence.

### CT mode

In `CT` mode:

```text
Reference C matches query C or T
```

Examples:

```text
Reference C, query C  -> match
Reference C, query T  -> match
Reference T, query C  -> mismatch
```

This directionality is important. The rule does not treat `C` and `T` as interchangeable in both directions.

### GA mode

In `GA` mode:

```text
Reference G matches query G or A
```

Examples:

```text
Reference G, query G  -> match
Reference G, query A  -> match
Reference A, query G  -> mismatch
```

### EITHER mode

In `EITHER` mode, the program performs both CT-aware and GA-aware alignments and reports the higher-scoring result.

If the two modes produce equal scores, CT mode is selected.

## Why a custom implementation?

Established tools remain the preferred choice for their intended applications.

### Bismark

Bismark is designed for genome-wide alignment and methylation analysis of bisulfite-sequencing reads. It supports workflows involving:

* reference-genome preparation;
* bisulfite-converted genome indexes;
* FASTQ read mapping;
* alignment backends such as Bowtie 2;
* SAM or BAM output;
* methylation extraction;
* strand and conversion-state interpretation;
* experiment-level quality-control reports.

Those capabilities are essential for whole-genome bisulfite-sequencing analysis but are unnecessary when the task is simply to align two already selected sequences.

Using Bismark for this restricted pairwise comparison would require substantially more setup and output processing without changing the underlying pairwise-alignment question.

### Parasail

Parasail is a high-performance pairwise sequence-alignment library with SIMD-accelerated implementations of several alignment algorithms.

Parasail is appropriate when alignment speed, large-scale processing, or integration with an optimized alignment library is required. For the limited purpose addressed by `nw_bs`, using Parasail would still require additional code to:

* define a bisulfite-specific substitution matrix;
* select the desired alignment routine;
* configure gap penalties;
* extract traceback information;
* interpret library result objects;
* format the alignment;
* generate the exact text fields needed by downstream scripts.

`nw_bs` instead keeps the complete scoring, matching, traceback, and output behavior in one short source file.

The reason for using `nw_bs` is therefore not that Bismark or Parasail are inadequate. They address broader or more performance-intensive problems. For this narrowly defined task, they add complexity that is not required.

## Advantages of nw_bs

The principal advantages are:

* **Restricted and explicit behavior**
  The program performs one clearly defined operation.

* **Directed bisulfite matching**
  CT and GA conversion rules are implemented directly in the substitution function.

* **Deterministic output**
  Traceback ties are resolved in a fixed order:

  ```text
  diagonal > up > left
  ```

* **Minimal dependencies**
  Only a C compiler and the standard C library are required.

* **Easy inspection**
  The complete algorithm is contained in one source file.

* **Easy customization**
  Matching rules, scoring, tie-breaking, and output fields can be modified directly.

* **Easy parsing**
  The default output is compact and line-oriented.

* **Reproducibility**
  The exact alignment behavior does not depend on external indexes, wrapper scripts, or library-specific result structures.

## Compilation

Compile with GCC:

```bash
gcc -O2 -Wall -Wextra -o nw_bs nw_bs.c
```

A stricter build may be performed with:

```bash
gcc \
    -std=c11 \
    -O2 \
    -Wall \
    -Wextra \
    -Wpedantic \
    -o nw_bs \
    nw_bs.c
```

No external libraries are required.

## Usage

```text
nw_bs --mode CT|GA|EITHER [options] REFSEQ READSEQ
```

Options:

```text
--mode CT|GA|EITHER
    Select the bisulfite-aware alignment mode.

-m INT
    Match score.

-x INT
    Mismatch score.

-g INT
    Linear gap score.
```

Default scores:

```text
match     =  1
mismatch  = -1
gap       = -1
```

## Examples

### CT-aware alignment

```bash
./nw_bs --mode CT ACGCCG ACTTCG
```

Reference cytosines aligned to query thymines are treated as compatible matches.

### GA-aware alignment

```bash
./nw_bs --mode GA GGGCCC AAACCC
```

Reference guanines aligned to query adenines are treated as compatible matches.

### Select the better conversion mode

```bash
./nw_bs --mode EITHER GGGCCC AAACCC
```

The program performs both CT-aware and GA-aware alignments and reports the higher-scoring result.

### Custom scoring

```bash
./nw_bs \
    --mode EITHER \
    -m 2 \
    -x -1 \
    -g -2 \
    REFSEQ \
    READSEQ
```

## Output

The output has five lines:

```text
mode:  CT
score: 12
ACGCCG
||||||
ACTTCG
```

The fields are:

1. selected bisulfite mode;
2. alignment score;
3. aligned reference sequence;
4. match indicator line;
5. aligned query sequence.

The compact output was chosen to support both direct inspection and simple downstream parsing.

For example:

```bash
mode=$(./nw_bs --mode EITHER REFSEQ READSEQ | awk '/^mode:/ {print $2}')
```

or:

```bash
score=$(./nw_bs --mode EITHER REFSEQ READSEQ | awk '/^score:/ {print $2}')
```

The output format can be readily modified to produce:

* tab-separated fields;
* CSV output;
* JSON output;
* CIGAR-like strings;
* mismatch counts;
* bisulfite-conversion counts;
* exact-match counts;
* gap counts;
* percent identity;
* machine-readable one-line summaries.

## Algorithm

The implementation uses the standard Needleman–Wunsch global-alignment recurrence with a linear gap penalty.

For each matrix cell, three possible transitions are evaluated:

```text
diagonal: previous diagonal score + substitution score
up:       score above + gap penalty
left:     score to the left + gap penalty
```

The highest-scoring transition is retained.

When multiple transitions have the same score, traceback uses the following fixed priority:

```text
DIAG > UP > LEFT
```

This produces one deterministic optimal alignment.

## Global-alignment behavior

`nw_bs` performs global alignment.

This means the complete reference sequence and complete query sequence are included in the alignment. Terminal gaps are penalized in the same manner as internal gaps.

This behavior is appropriate when the two input sequences are expected to represent corresponding full-length regions.

It may not be appropriate when:

* the reference is a long genomic interval containing a shorter read;
* only a subsequence is expected to align;
* unaligned terminal sequence should not be penalized;
* local or semiglobal alignment is required.

For those applications, the initialization and traceback rules should be modified or a more suitable alignment program should be used.

## Gap model

The current implementation uses a linear gap penalty.

A gap of length (k) receives:

```text
k × gap_score
```

The program does not currently implement separate gap-opening and gap-extension penalties.

For analyses in which biologically realistic insertion and deletion modeling is important, an affine-gap implementation such as the Gotoh algorithm may be preferable.

## Treatment of N

By default, `N` in either sequence is treated as matching any base.

This behavior is defined in `bs_match()` and can be changed easily.

For example, an application may instead choose to treat `N` as:

* a mismatch;
* a neutral score of zero;
* a match only when aligned to another `N`.

Users should select behavior appropriate for their data and report any modification.

## Computational complexity

For reference length (n) and query length (m), the implementation requires:

```text
Time:   O(nm)
Memory: O(nm)
```

The full dynamic-programming and traceback matrices are retained in memory.

This is appropriate for short predefined sequence pairs. It is not intended for chromosome-scale alignment or high-throughput mapping of millions of reads.

## Customization

Most behavior can be changed in a small number of functions.

### Bisulfite matching

Modify:

```c
static inline int bs_match(...)
```

This function defines which reference/query base combinations are considered compatible.

### Substitution scoring

Modify:

```c
static inline int sub_score(...)
```

This function assigns match and mismatch scores.

### Traceback tie-breaking

Modify the direction-selection block in the dynamic-programming loop.

The current priority is:

```text
DIAG > UP > LEFT
```

### Output formatting

Modify the final `printf()` statements in `main()`.

Because the output logic is independent of the dynamic-programming calculation, it can be adapted for a specific downstream workflow without changing the alignment algorithm.

## Validation suggestions

Useful tests include the following.

### Exact match

```bash
./nw_bs --mode CT ACGT ACGT
```

### CT-compatible conversion

```bash
./nw_bs --mode CT CCCC TTTT
```

### CT directionality

```bash
./nw_bs --mode CT TTTT CCCC
```

The reverse pairing should not be treated as a bisulfite-compatible match.

### GA-compatible conversion

```bash
./nw_bs --mode GA GGGG AAAA
```

### EITHER mode selection

```bash
./nw_bs --mode EITHER GGGG AAAA
```

GA mode should be selected.

### Gap placement

```bash
./nw_bs --mode CT ACGTT ACGT
```

### Lowercase input

```bash
./nw_bs --mode CT acgccg acttcg
```

Input comparison is case-insensitive.

## Reproducibility

For published analyses, report:

* repository URL;
* release or commit identifier;
* compiler and compiler version;
* compilation command;
* selected mode;
* match score;
* mismatch score;
* gap score;
* treatment of `N`;
* exact command line used.

A versioned release corresponding to a manuscript should preferably be archived with a permanent identifier through a service such as Zenodo.

## Citation

When using `nw_bs`, cite the associated publication or repository release:

```text
[Citation to be added]
```

## License

Add the selected software license here.

A permissive license such as MIT, BSD-2-Clause, or BSD-3-Clause is often suitable for a small research utility, subject to institutional requirements.

## Author

```text
[Author name]
[Institution]
[Contact information]
```

## Disclaimer

`nw_bs` is research software intended for a narrowly defined pairwise-alignment task. Users are responsible for validating the scoring model, conversion rules, and output for their application. It should not be used as a substitute for established genome-scale bisulfite-sequencing alignment and methylation-analysis pipelines.
