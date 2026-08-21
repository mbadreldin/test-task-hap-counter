# Counts TSV Output Format

Each row in `counts.tsv` represents one eligible biallelic SNV from the input
VCF. Positions are one-based. Reads without a valid `HP=1` or `HP=2` tag do
not contribute to the haplotype counts.

| Column | Meaning |
| --- | --- |
| `chrom` | Chromosome or contig name from the VCF, such as `chr16`. |
| `pos` | One-based genomic position of the SNV. |
| `h1_REF` | Number of `HP=1` alignments whose observed base matches the VCF reference allele. |
| `h1_ALT` | Number of `HP=1` alignments whose observed base matches the VCF alternate allele. |
| `h2_REF` | Number of `HP=2` alignments whose observed base matches the VCF reference allele. |
| `h2_ALT` | Number of `HP=2` alignments whose observed base matches the VCF alternate allele. |
| `h1_call` | Inferred haplotype 1 allele: `REF`, `ALT`, or `NO_CALL`. |
| `h1_discrepant` | Minority-allele observations for haplotype 1: `min(h1_REF, h1_ALT)`. |
| `h1_discrepancy_fraction` | Haplotype 1 minority fraction: `h1_discrepant / (h1_REF + h1_ALT)`. |
| `h2_call` | Inferred haplotype 2 allele: `REF`, `ALT`, or `NO_CALL`. |
| `h2_discrepant` | Minority-allele observations for haplotype 2: `min(h2_REF, h2_ALT)`. |
| `h2_discrepancy_fraction` | Haplotype 2 minority fraction: `h2_discrepant / (h2_REF + h2_ALT)`. |

The first six columns are always present. The call and discrepancy columns are
added when `--genotype` or `--histogram` is requested. `NO_CALL` means the
evidence did not meet the configured depth or winning-allele-fraction
threshold. Discrepancy values remain available for covered `NO_CALL` cells;
they are blank only when the informative depth is zero.

## How Haplotype Calls Are Inferred

Each haplotype is evaluated independently. For haplotype 1, the informative
depth is:

```text
h1_depth = h1_REF + h1_ALT
```

Only observations matching the VCF REF or ALT allele contribute to this depth.
The winning allele is whichever has the larger count, and its fraction is:

```text
winning_fraction = max(h1_REF, h1_ALT) / h1_depth
```

The program reports `REF` or `ALT` only when all of these conditions hold:

1. Depth is at least `--min-depth` (default: `5`).
2. REF and ALT counts are not tied.
3. The winning fraction is at least `--min-call-fraction` (default: `0.70`).

Otherwise, it reports `NO_CALL`. For example, counts of 8 REF and 2 ALT have a
depth of 10 and a winning fraction of 0.80, producing `REF` with the defaults.
Counts of 6 REF and 4 ALT have sufficient depth but a winning fraction of only
0.60, producing `NO_CALL`. The same calculation is applied to haplotype 2.

Depth measures the amount of informative evidence, while the winning fraction
measures how consistently that evidence supports one allele. Configure them
independently with `--min-depth` and `--min-call-fraction`.

## Interpreting the Histogram

The optional histogram plots one discrepancy fraction for every covered
haplotype/SNV cell, including covered cells that resulted in `NO_CALL`:

```text
discrepancy_fraction = min(REF, ALT) / (REF + ALT)
```

The x-axis ranges from 0 to 0.5. Values near 0 mean nearly all observations
agree on one allele; values near 0.5 mean REF and ALT observations are evenly
mixed. The y-axis is the number of covered haplotype/SNV cells in each range.
A distribution concentrated near 0 therefore indicates consistent
haplotype-specific allele support. A substantial tail toward 0.5 indicates more
mixed support, which may reflect sequencing errors, alignment ambiguity,
incorrect haplotags, or genuinely difficult loci.

The histogram summarizes consistency, not coverage: a low-depth cell and a
high-depth cell each contribute one value. It also does not distinguish called
cells from `NO_CALL` cells, so use `counts.tsv` and `summary.txt` when assessing
callability or evidence depth.

## Example

```text
chrom  pos       h1_REF  h1_ALT  h2_REF  h2_ALT  h1_call  h1_discrepant  h1_discrepancy_fraction  h2_call  h2_discrepant  h2_discrepancy_fraction
chr16  28001381  11      4       3       12      REF      4               0.266667                 ALT      3               0.200000
```

At `chr16:28001381`, haplotype 1 has 11 REF and 4 ALT observations, so it is
called `REF`; its discrepancy fraction is `4 / 15 = 0.266667`. Haplotype 2 has
3 REF and 12 ALT observations, so it is called `ALT`; its discrepancy fraction
is `3 / 15 = 0.200000`.
