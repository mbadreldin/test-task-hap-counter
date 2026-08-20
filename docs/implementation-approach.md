# Implementation Approach

## Decision

The counter uses a variant-by-variant pileup. It reads eligible SNVs from the VCF in file order and requests a one-base pileup from the indexed BAM for each variant. Pysam/htslib performs the CIGAR-aware mapping from reference position to read base.

This approach is compact, easy to explain, and gives each VCF record a direct path to one TSV row. Filtering primary alignments, checking `HP`, and classifying REF or ALT observations remain visible in application code rather than being hidden behind a larger traversal abstraction.

Only distinct A/C/G/T biallelic SNVs are eligible. VCF contigs are checked against the BAM header before counting so chromosome naming mistakes fail clearly instead of producing plausible-looking zero counts. Missing and invalid HP tags do not contribute to required counts, but their informative REF/ALT observations are retained as run diagnostics.

## Performance Trade-off

The supplied data contains 939 relevant SNVs and 3,269 alignments. A read-only benchmark on the development machine counted 22,107 informative votes in approximately 11 seconds. That is adequate for the task and keeps implementation risk low within the suggested two-hour scope.

A chromosome-streamed prototype produced the same votes in approximately 0.8 seconds. It avoids repeated BAM seeks and decoding of blocks containing long reads that overlap multiple variants. Its cost is additional position-indexing and result-retention logic, and it traverses pileup columns that may not contain variants.

The timings are illustrative, not portable benchmarks. For whole-genome VCFs or workloads where repeated queries dominate runtime, the implementation should move to chromosome-streamed pileup or a read/variant interval join. The counting and result models are kept separate so that traversal can be replaced without changing output or genotype analysis.

## Configurable Genotyping

“Enough haplotagged alignments” is defined independently for each haplotype as its informative depth: `REF + ALT` votes. `--min-depth` controls the minimum and defaults to 5.

Coverage alone does not show whether a call is decisive, so `--min-call-fraction` separately controls the required winning-allele fraction and defaults to 0.70. A call is `NO_CALL` when depth is insufficient, counts are tied, or the winning fraction is below the threshold. Keeping both values configurable makes the policy explicit and allows adjustment for sequencing depth and error rate without changing code.

Discrepancy is the minority of REF and ALT votes divided by informative depth. The summary reports it for all covered haplotype/SNV cells and separately for cells that pass calling thresholds. The histogram uses all covered cells and the mathematically possible range from 0 to 0.5.
