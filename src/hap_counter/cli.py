"""Command-line interface for hap-counter."""

import argparse
import sys
from collections.abc import Sequence
from contextlib import ExitStack
from pathlib import Path

from hap_counter.analysis import summarize_analysis
from hap_counter.counter import ReadFilterConfig, count_haplotype_support
from hap_counter.output import (
    plot_discrepancy_histogram,
    write_analysis_summary,
    write_support_tsv,
)
from hap_counter.vcf import load_biallelic_snvs


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hap-counter",
        description="Count haplotype-specific REF and ALT support at biallelic SNVs.",
    )
    parser.add_argument("bam", type=Path, help="coordinate-sorted, indexed BAM")
    parser.add_argument("vcf", type=Path, help="plain or compressed VCF")
    parser.add_argument("-o", "--output", required=True, help="output TSV path or '-'")
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--min-baseq", type=int, default=0)
    parser.add_argument("--genotype", action="store_true", help="append allele calls")
    parser.add_argument(
        "--min-depth",
        type=int,
        default=5,
        help="minimum REF+ALT votes per haplotype call (default: 5)",
    )
    parser.add_argument(
        "--min-call-fraction",
        type=float,
        default=0.70,
        help="minimum winning-allele fraction (default: 0.70)",
    )
    parser.add_argument("--histogram", type=Path, help="optional discrepancy PNG")
    return parser


def main(command_line_arguments: Sequence[str] | None = None) -> int:
    parser = build_argument_parser()
    arguments = parser.parse_args(command_line_arguments)
    _validate_cli_arguments(parser, arguments)

    try:
        snvs = load_biallelic_snvs(arguments.vcf)
        support_counts = count_haplotype_support(
            arguments.bam,
            snvs,
            ReadFilterConfig(
                minimum_mapping_quality=arguments.min_mapq,
                minimum_base_quality=arguments.min_baseq,
            ),
        )
        include_analysis = arguments.genotype or arguments.histogram is not None
        with ExitStack() as output_stack:
            output_stream = sys.stdout
            if arguments.output != "-":
                output_stream = output_stack.enter_context(
                    open(arguments.output, "w", encoding="utf-8", newline="")
                )
            calls_by_snv = write_support_tsv(
                support_counts,
                output_stream,
                include_analysis=include_analysis,
                minimum_depth=arguments.min_depth,
                minimum_call_fraction=arguments.min_call_fraction,
            )

        if arguments.histogram is not None:
            plot_discrepancy_histogram(calls_by_snv, arguments.histogram)

        if include_analysis:
            summary = summarize_analysis(support_counts, calls_by_snv)
            write_analysis_summary(
                summary,
                sys.stderr,
                eligible_snv_count=len(snvs),
                tsv_output=arguments.output,
                histogram_output=(
                    str(arguments.histogram)
                    if arguments.histogram is not None
                    else None
                ),
            )
        else:
            print(f"eligible_snvs={len(snvs)}", file=sys.stderr)
            print(f"wrote_tsv={arguments.output}", file=sys.stderr)
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(2, f"hap-counter: error: {error}\n")
    return 0


def _validate_cli_arguments(
    parser: argparse.ArgumentParser, arguments: argparse.Namespace
) -> None:
    if arguments.min_mapq < 0 or arguments.min_baseq < 0:
        parser.error("quality thresholds must be non-negative")
    if arguments.min_depth < 1:
        parser.error("--min-depth must be at least 1")
    if not 0.5 <= arguments.min_call_fraction <= 1.0:
        parser.error("--min-call-fraction must be between 0.5 and 1.0")
