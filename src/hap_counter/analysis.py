"""Optional genotype and discrepancy analysis."""

from collections.abc import Iterable
from dataclasses import dataclass

from hap_counter.models import HaplotypeSupportCounts


@dataclass(frozen=True, slots=True)
class HaplotypeCall:
    """An inferred allele and the observations behind it."""

    allele: str
    informative_depth: int
    winning_allele_fraction: float | None
    discrepant_votes: int | None
    discrepancy_fraction: float | None

    @property
    def is_called(self) -> bool:
        return self.allele != "NO_CALL"


@dataclass(frozen=True, slots=True)
class DiscrepancyMetrics:
    """Aggregated discrepancy statistics for haplotype/SNV cells."""

    cell_count: int
    informative_votes: int
    discrepant_votes: int

    @property
    def discrepancy_fraction(self) -> float:
        if self.informative_votes == 0:
            return 0.0
        return self.discrepant_votes / self.informative_votes


@dataclass(frozen=True, slots=True)
class AnalysisSummary:
    """Run-level coverage, discrepancy, and untagged-read diagnostics."""

    covered: DiscrepancyMetrics
    called: DiscrepancyMetrics
    untagged_ref_votes: int
    untagged_alt_votes: int

    @property
    def untagged_informative_votes(self) -> int:
        return self.untagged_ref_votes + self.untagged_alt_votes


PerSnvCalls = tuple[HaplotypeCall, HaplotypeCall]


def infer_haplotype_call(
    ref_count: int,
    alt_count: int,
    *,
    minimum_depth: int = 5,
    minimum_call_fraction: float = 0.70,
) -> HaplotypeCall:
    """Infer REF or ALT using configurable depth and allele-fraction thresholds."""
    if minimum_depth < 1:
        raise ValueError("minimum_depth must be at least 1")
    if not 0.5 <= minimum_call_fraction <= 1.0:
        raise ValueError("minimum_call_fraction must be between 0.5 and 1.0")
    if ref_count < 0 or alt_count < 0:
        raise ValueError("allele counts must be non-negative")

    informative_depth = ref_count + alt_count
    if informative_depth == 0:
        return HaplotypeCall("NO_CALL", 0, None, None, None)

    winning_count = max(ref_count, alt_count)
    discrepant_votes = min(ref_count, alt_count)
    winning_fraction = winning_count / informative_depth
    discrepancy_fraction = discrepant_votes / informative_depth
    enough_evidence = (
        informative_depth >= minimum_depth
        and ref_count != alt_count
        and winning_fraction >= minimum_call_fraction
    )
    allele = "REF" if ref_count > alt_count else "ALT"
    if not enough_evidence:
        allele = "NO_CALL"

    return HaplotypeCall(
        allele=allele,
        informative_depth=informative_depth,
        winning_allele_fraction=winning_fraction,
        discrepant_votes=discrepant_votes,
        discrepancy_fraction=discrepancy_fraction,
    )


def summarize_analysis(
    support_counts: Iterable[HaplotypeSupportCounts],
    calls_by_snv: Iterable[PerSnvCalls],
) -> AnalysisSummary:
    """Aggregate covered-cell, called-cell, and untagged-read statistics."""
    covered_cells = covered_votes = covered_discrepant = 0
    called_cells = called_votes = called_discrepant = 0
    untagged_ref_votes = untagged_alt_votes = 0

    for support, haplotype_calls in zip(
        support_counts, calls_by_snv, strict=True
    ):
        untagged_ref_votes += support.untagged_ref
        untagged_alt_votes += support.untagged_alt
        for haplotype_call in haplotype_calls:
            if haplotype_call.informative_depth > 0:
                covered_cells += 1
                covered_votes += haplotype_call.informative_depth
                covered_discrepant += haplotype_call.discrepant_votes or 0
            if haplotype_call.is_called:
                called_cells += 1
                called_votes += haplotype_call.informative_depth
                called_discrepant += haplotype_call.discrepant_votes or 0

    return AnalysisSummary(
        covered=DiscrepancyMetrics(
            covered_cells, covered_votes, covered_discrepant
        ),
        called=DiscrepancyMetrics(called_cells, called_votes, called_discrepant),
        untagged_ref_votes=untagged_ref_votes,
        untagged_alt_votes=untagged_alt_votes,
    )
