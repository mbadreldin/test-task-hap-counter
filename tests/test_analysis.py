import pytest

from hap_counter.analysis import infer_haplotype_call, summarize_analysis
from hap_counter.models import BiallelicSnv, HaplotypeSupportCounts


def test_call_requires_minimum_depth() -> None:
    call = infer_haplotype_call(
        3, 1, minimum_depth=5, minimum_call_fraction=0.7
    )

    assert call.allele == "NO_CALL"
    assert call.informative_depth == 4
    assert call.discrepancy_fraction == pytest.approx(0.25)


def test_call_requires_configured_winning_fraction() -> None:
    call = infer_haplotype_call(
        6, 4, minimum_depth=5, minimum_call_fraction=0.7
    )

    assert call.allele == "NO_CALL"
    assert call.winning_allele_fraction == pytest.approx(0.6)
    assert call.discrepant_votes == 4


def test_call_reports_discrepant_votes() -> None:
    call = infer_haplotype_call(
        2, 8, minimum_depth=5, minimum_call_fraction=0.7
    )

    assert call.allele == "ALT"
    assert call.discrepant_votes == 2
    assert call.discrepancy_fraction == pytest.approx(0.2)


def test_tie_is_not_called_even_at_half_fraction() -> None:
    call = infer_haplotype_call(
        3, 3, minimum_depth=5, minimum_call_fraction=0.5
    )

    assert call.allele == "NO_CALL"
    assert call.discrepancy_fraction == pytest.approx(0.5)


@pytest.mark.parametrize(
    ("minimum_depth", "minimum_fraction"),
    [(0, 0.7), (5, 0.49), (5, 1.01)],
)
def test_invalid_call_thresholds_are_rejected(
    minimum_depth: int, minimum_fraction: float
) -> None:
    with pytest.raises(ValueError):
        infer_haplotype_call(
            5,
            0,
            minimum_depth=minimum_depth,
            minimum_call_fraction=minimum_fraction,
        )


def test_summary_reports_covered_called_and_untagged_metrics() -> None:
    support = HaplotypeSupportCounts(
        snv=BiallelicSnv("chr1", 101, "C", "G"),
        haplotype_1_ref=8,
        haplotype_1_alt=2,
        haplotype_2_ref=1,
        haplotype_2_alt=1,
        untagged_ref=3,
        untagged_alt=4,
    )
    calls = (
        infer_haplotype_call(8, 2),
        infer_haplotype_call(1, 1),
    )

    summary = summarize_analysis([support], [calls])

    assert summary.covered.cell_count == 2
    assert summary.covered.informative_votes == 12
    assert summary.covered.discrepant_votes == 3
    assert summary.called.cell_count == 1
    assert summary.called.informative_votes == 10
    assert summary.called.discrepant_votes == 2
    assert summary.untagged_ref_votes == 3
    assert summary.untagged_alt_votes == 4
    assert summary.untagged_informative_votes == 7
