# How to Run

## Install

From the repository root, create an isolated environment and install the project. The `analysis` extra provides Matplotlib for histogram output.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[analysis]'
.venv/bin/python -c "import hap_counter, pysam, matplotlib"
.venv/bin/hap-counter --help
```

Using the executable inside `.venv` avoids accidentally running with packages installed in a different Python environment. It also prevents a local `hap_counter.py` file from shadowing the installed package.

The input BAM must be coordinate-sorted and have a matching `.bai` index. The VCF may be plain text or bgzip-compressed.

## Counts, Genotypes, and Histogram

Run the complete analysis on the supplied data with:

```bash
mkdir -p output
MPLCONFIGDIR=output/.matplotlib \
XDG_CACHE_HOME=output/.cache \
.venv/bin/hap-counter \
  test_data/giab_2023.05.hg002.haplotagged.chr16_28000000_29000000.processed.30x.bam \
  test_data/giab_2023.05.hg002.wf_snp.chr16_28000000_29000000.vcf.gz \
  --output output/counts.tsv \
  --genotype \
  --min-depth 5 \
  --min-call-fraction 0.70 \
  --histogram output/discrepancies.png \
  2> output/summary.txt
```

This creates:

- `output/counts.tsv`: mandatory HP1/HP2 REF and ALT counts, inferred haplotype calls, and per-call discrepancy values.
- `output/discrepancies.png`: histogram of discrepancy fractions for all covered haplotype/SNV cells.
- `output/summary.txt`: eligible-SNV, untagged-read, covered-cell, and called-cell statistics.

The two cache variables keep Matplotlib and font-cache files under the ignored `output/` directory, avoiding environment warnings in the summary.

`--min-depth` defines how many informative `REF + ALT` observations a haplotype needs. `--min-call-fraction` defines the minimum fraction supporting its winning allele. Adjust both for the desired calling stringency. Untagged votes are diagnostics only and never enter HP1/HP2 calls. `output/summary.txt` contains application statistics only after the command starts successfully; if output files are missing, run the help preflight above without redirection.

## Mandatory Counts Only

Run the mandatory task on the supplied data with:

```bash
mkdir -p output
.venv/bin/hap-counter \
  test_data/giab_2023.05.hg002.haplotagged.chr16_28000000_29000000.processed.30x.bam \
  test_data/giab_2023.05.hg002.wf_snp.chr16_28000000_29000000.vcf.gz \
  --output output/counts.tsv
```

This produces the six required TSV columns without inferred genotype calls, discrepancy analysis, or a histogram. Use `--output -` to write TSV data to standard output. Run `.venv/bin/hap-counter --help` for all filtering options.
