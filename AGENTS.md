# Repository Guidelines

## Project Structure & Module Organization

`README.md` is the authoritative task specification. Python code belongs in `src/hap_counter/`, automated tests belong in `tests/`, and design decisions are documented in `docs/`. Reference genomic fixtures live in `test_data/`: the coordinate-sorted BAM has a matching `.bai`, and the phased VCF is available as compressed `.vcf.gz` data with a `.tbi` index.

Keep allele-counting logic separate from BAM/VCF parsing so each part can be tested independently. Do not edit reference fixtures in place.

## Build, Test, and Development Commands

Create a local environment and install the package with development tools:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[dev]'
.venv/bin/python -m pytest
.venv/bin/python -m ruff check src tests
```

The install command adds the package in editable mode. Pytest runs the suite, while Ruff checks style without rewriting files. In PyCharm, open the repository root, choose `.venv/bin/python` as the project interpreter, and mark `src` as a Sources Root if it is not detected automatically.

Run the program with `.venv/bin/hap-counter INPUT.bam INPUT.vcf.gz -o counts.tsv`. Add `--genotype` for configurable haplotype calls or `--histogram discrepancies.png` for optional analysis.

## Coding Style & Naming Conventions

Use four-space indentation and follow PEP 8. Name functions and variables with `snake_case`, classes with `PascalCase`, and constants with `UPPER_SNAKE_CASE`. Add type hints to public functions and concise docstrings where behavior is not obvious. Keep functions small, especially around coordinate conversion, CIGAR traversal, and allele classification.

## Testing Guidelines

Cover REF and ALT observations on both haplotypes, missing or invalid `HP` tags, secondary and supplementary alignments, uncovered variants, and CIGAR operations that alter read/reference coordinates. Test compressed and uncompressed VCF inputs. Use descriptive names such as `test_secondary_alignment_is_ignored` and small synthetic fixtures where practical.

## Commit & Pull Request Guidelines

Use short, imperative commit subjects, such as `Handle missing HP tags`. Pull requests should explain the behavior and assumptions, link relevant issues, list verification commands, and include representative TSV output when output behavior changes.

Do not commit caches, generated output, unnecessary large genomic datasets, or sensitive sample data.
