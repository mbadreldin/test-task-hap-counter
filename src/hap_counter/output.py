"""TSV, summary, and histogram output."""

import csv
from collections.abc import Iterable
from pathlib import Path
from typing import TextIO

from hap_counter.analysis import (
    AnalysisSummary,
    DiscrepancyMetrics,
    HaplotypeCall,
    PerSnvCalls,
    infer_haplotype_call,
)
from hap_counter.models import HaplotypeSupportCounts


def write_support_tsv(
    support_counts: Iterable[HaplotypeSupportCounts],
    output_stream: TextIO,
    *,
    include_analysis: bool = False,
    minimum_depth: int = 5,
    minimum_call_fraction: float = 0.70,
) -> list[PerSnvCalls]:
    """Write required counts and, optionally, call/discrepancy columns."""
    header = ["chrom", "pos", "h1_REF", "h1_ALT", "h2_REF", "h2_ALT"]
    if include_analysis:
        for haplotype in (1, 2):
            header.extend(
                [
                    f"h{haplotype}_call",
                    f"h{haplotype}_discrepant",
                    f"h{haplotype}_discrepancy_fraction",
                ]
            )

    writer = csv.writer(output_stream, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    calls_by_snv: list[PerSnvCalls] = []
    for support in support_counts:
        snv = support.snv
        row: list[str | int] = [
            snv.chromosome,
            snv.position,
            support.haplotype_1_ref,
            support.haplotype_1_alt,
            support.haplotype_2_ref,
            support.haplotype_2_alt,
        ]
        if include_analysis:
            haplotype_calls = (
                infer_haplotype_call(
                    *support.counts_for_haplotype(1),
                    minimum_depth=minimum_depth,
                    minimum_call_fraction=minimum_call_fraction,
                ),
                infer_haplotype_call(
                    *support.counts_for_haplotype(2),
                    minimum_depth=minimum_depth,
                    minimum_call_fraction=minimum_call_fraction,
                ),
            )
            calls_by_snv.append(haplotype_calls)
            for haplotype_call in haplotype_calls:
                row.extend(_format_haplotype_call(haplotype_call))
        writer.writerow(row)
    return calls_by_snv


def _format_haplotype_call(haplotype_call: HaplotypeCall) -> list[str | int]:
    return [
        haplotype_call.allele,
        ""
        if haplotype_call.discrepant_votes is None
        else haplotype_call.discrepant_votes,
        ""
        if haplotype_call.discrepancy_fraction is None
        else f"{haplotype_call.discrepancy_fraction:.6f}",
    ]


def write_analysis_summary(
    summary: AnalysisSummary,
    output_stream: TextIO,
    *,
    eligible_snv_count: int,
    tsv_output: str,
    histogram_output: str | None,
) -> None:
    """Write concise run statistics as key-value lines."""
    print(f"eligible_snvs={eligible_snv_count}", file=output_stream)
    print(
        f"untagged_ref_votes={summary.untagged_ref_votes}\t"
        f"untagged_alt_votes={summary.untagged_alt_votes}\t"
        f"untagged_informative_votes={summary.untagged_informative_votes}",
        file=output_stream,
    )
    _write_discrepancy_metrics("covered", summary.covered, output_stream)
    _write_discrepancy_metrics("called", summary.called, output_stream)
    print(f"wrote_tsv={tsv_output}", file=output_stream)
    if histogram_output is not None:
        print(f"wrote_histogram={histogram_output}", file=output_stream)


def _write_discrepancy_metrics(
    prefix: str, metrics: DiscrepancyMetrics, output_stream: TextIO
) -> None:
    print(
        f"{prefix}_cells={metrics.cell_count}\t"
        f"{prefix}_votes={metrics.informative_votes}\t"
        f"{prefix}_discrepant_votes={metrics.discrepant_votes}\t"
        f"{prefix}_discrepancy_fraction={metrics.discrepancy_fraction:.6f}",
        file=output_stream,
    )


def plot_discrepancy_histogram(
    calls_by_snv: Iterable[PerSnvCalls], output_path: str | Path
) -> None:
    """Plot discrepancy fractions for all covered haplotype/SNV cells."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError(
            "histogram output requires the 'analysis' extra: "
            "python -m pip install -e '.[analysis]'"
        ) from error

    discrepancy_fractions = [
        call.discrepancy_fraction
        for haplotype_calls in calls_by_snv
        for call in haplotype_calls
        if call.discrepancy_fraction is not None
    ]
    figure, axis = plt.subplots(figsize=(8, 5))
    axis.hist(
        discrepancy_fractions,
        bins=20,
        range=(0.0, 0.5),
        edgecolor="black",
    )
    axis.set_xlim(0.0, 0.5)
    axis.set_xlabel("Discrepant vote fraction (minority allele / informative depth)")
    axis.set_ylabel("Covered haplotype/SNV cells")
    axis.set_title(f"Haplotype vote discrepancies (n={len(discrepancy_fractions)})")
    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)
