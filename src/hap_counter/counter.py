"""CIGAR-aware, variant-by-variant allele counting."""

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import pysam

from hap_counter.models import BiallelicSnv, HaplotypeSupportCounts


@dataclass(frozen=True, slots=True)
class ReadFilterConfig:
    """Quality thresholds applied before an alignment contributes a vote."""

    minimum_mapping_quality: int = 0
    minimum_base_quality: int = 0

    def __post_init__(self) -> None:
        if self.minimum_mapping_quality < 0 or self.minimum_base_quality < 0:
            raise ValueError("quality thresholds must be non-negative")


def count_haplotype_support(
    bam_path: str | Path,
    snvs: Iterable[BiallelicSnv],
    read_filters: ReadFilterConfig | None = None,
) -> list[HaplotypeSupportCounts]:
    """Count haplotype-specific REF and ALT observations for each SNV."""
    filters = read_filters or ReadFilterConfig()
    snv_list = list(snvs)

    with pysam.AlignmentFile(str(bam_path), "rb") as alignments:
        _validate_bam_contigs(snv_list, set(alignments.references))
        return [
            _count_support_for_snv(alignments, snv, filters) for snv in snv_list
        ]


def _validate_bam_contigs(
    snvs: Iterable[BiallelicSnv], bam_contigs: set[str]
) -> None:
    missing_contigs = sorted({snv.chromosome for snv in snvs} - bam_contigs)
    if missing_contigs:
        formatted_contigs = ", ".join(missing_contigs)
        raise ValueError(f"VCF contigs are absent from the BAM header: {formatted_contigs}")


def _count_support_for_snv(
    alignments: pysam.AlignmentFile,
    snv: BiallelicSnv,
    read_filters: ReadFilterConfig,
) -> HaplotypeSupportCounts:
    support_counts = HaplotypeSupportCounts(snv=snv)
    # VCF positions are one-based; pysam pileup coordinates are zero-based.
    reference_position = snv.position - 1
    pileup_columns = alignments.pileup(
        snv.chromosome,
        reference_position,
        reference_position + 1,
        truncate=True,
        # Apply only the task's explicit primary-alignment rules below.
        stepper="nofilter",
        min_base_quality=0,
        min_mapping_quality=0,
        ignore_overlaps=False,
        ignore_orphans=False,
        max_depth=1_000_000,
    )

    for pileup_column in pileup_columns:
        if pileup_column.reference_pos != reference_position:
            continue
        for pileup_read in pileup_column.pileups:
            alignment = pileup_read.alignment
            if not _alignment_can_vote(alignment, pileup_read, read_filters):
                continue

            query_position = pileup_read.query_position
            assert query_position is not None  # Guarded by _alignment_can_vote.
            query_sequence = alignment.query_sequence
            if query_sequence is None:
                continue

            observed_allele = _classify_observed_base(
                query_sequence[query_position].upper(), snv
            )
            if observed_allele is None:
                continue

            haplotype = _get_haplotype(alignment)
            if haplotype is None:
                support_counts.add_untagged_vote(observed_allele)
            else:
                support_counts.add_haplotype_vote(haplotype, observed_allele)

    return support_counts


def _alignment_can_vote(
    alignment: pysam.AlignedSegment,
    pileup_read: pysam.PileupRead,
    read_filters: ReadFilterConfig,
) -> bool:
    if (
        alignment.is_unmapped
        or alignment.is_secondary
        or alignment.is_supplementary
        or alignment.mapping_quality < read_filters.minimum_mapping_quality
        or pileup_read.is_del
        or pileup_read.is_refskip
        or pileup_read.query_position is None
    ):
        return False

    if read_filters.minimum_base_quality == 0:
        return True
    base_qualities = alignment.query_qualities
    return (
        base_qualities is not None
        and base_qualities[pileup_read.query_position]
        >= read_filters.minimum_base_quality
    )


def _classify_observed_base(base: str, snv: BiallelicSnv) -> str | None:
    if base == snv.reference_allele:
        return "REF"
    if base == snv.alternate_allele:
        return "ALT"
    # Deletions, skips, and third alleles support neither SNV allele.
    return None


def _get_haplotype(alignment: pysam.AlignedSegment) -> int | None:
    try:
        haplotype = alignment.get_tag("HP")
    except KeyError:
        return None
    return haplotype if haplotype in (1, 2) else None
