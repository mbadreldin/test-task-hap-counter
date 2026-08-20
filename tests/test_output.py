from io import StringIO
from pathlib import Path

from matplotlib.axes import Axes

from hap_counter.analysis import infer_haplotype_call, summarize_analysis
from hap_counter.counter import count_haplotype_support
from hap_counter.models import BiallelicSnv, HaplotypeSupportCounts
from hap_counter.output import (
    plot_discrepancy_histogram,
    write_analysis_summary,
    write_support_tsv,
)
from hap_counter.vcf import load_biallelic_snvs

FIXTURES = Path(__file__).parent / "fixtures"


def test_mandatory_output_has_six_columns() -> None:
    output_stream = StringIO()
    support = HaplotypeSupportCounts(
        snv=BiallelicSnv("chr1", 101, "C", "G"),
        haplotype_1_ref=2,
        haplotype_1_alt=1,
        haplotype_2_alt=3,
    )

    write_support_tsv([support], output_stream)

    assert output_stream.getvalue().splitlines() == [
        "chrom\tpos\th1_REF\th1_ALT\th2_REF\th2_ALT",
        "chr1\t101\t2\t1\t0\t3",
    ]


def test_optional_output_includes_calls_and_discrepancies() -> None:
    output_stream = StringIO()
    support = HaplotypeSupportCounts(
        snv=BiallelicSnv("chr1", 101, "C", "G"),
        haplotype_1_ref=8,
        haplotype_1_alt=2,
        haplotype_2_ref=1,
        haplotype_2_alt=9,
    )

    calls_by_snv = write_support_tsv(
        [support],
        output_stream,
        include_analysis=True,
        minimum_depth=5,
        minimum_call_fraction=0.7,
    )

    fields = output_stream.getvalue().splitlines()[1].split("\t")
    assert fields[6:] == ["REF", "2", "0.200000", "ALT", "1", "0.100000"]
    assert [call.allele for call in calls_by_snv[0]] == ["REF", "ALT"]


def test_optional_output_matches_reviewable_expected_file(edge_case_bam: Path) -> None:
    snvs = load_biallelic_snvs(FIXTURES / "variants.vcf")
    support_counts = count_haplotype_support(edge_case_bam, snvs)
    output_stream = StringIO()

    write_support_tsv(
        support_counts,
        output_stream,
        include_analysis=True,
        minimum_depth=1,
        minimum_call_fraction=0.7,
    )

    assert output_stream.getvalue() == (
        FIXTURES / "expected_analysis.tsv"
    ).read_text()


def test_histogram_uses_covered_cells_and_valid_fraction_range(
    tmp_path: Path, monkeypatch
) -> None:
    calls_by_snv = [
        (
            infer_haplotype_call(8, 2),
            infer_haplotype_call(1, 9),
        ),
        (
            infer_haplotype_call(1, 1),
            infer_haplotype_call(0, 0),
        ),
    ]
    captured_range = None
    original_hist = Axes.hist

    def capture_histogram_range(self, values, *args, **kwargs):
        nonlocal captured_range
        captured_range = kwargs.get("range")
        return original_hist(self, values, *args, **kwargs)

    monkeypatch.setattr(Axes, "hist", capture_histogram_range)
    output_path = tmp_path / "discrepancies.png"

    plot_discrepancy_histogram(calls_by_snv, output_path)

    assert captured_range == (0.0, 0.5)
    assert output_path.read_bytes().startswith(b"\x89PNG")


def test_summary_names_both_discrepancy_scopes_and_outputs() -> None:
    support = HaplotypeSupportCounts(
        snv=BiallelicSnv("chr1", 101, "C", "G"),
        haplotype_1_ref=8,
        haplotype_1_alt=2,
        haplotype_2_ref=1,
        haplotype_2_alt=1,
        untagged_alt=3,
    )
    calls_by_snv = [
        (
            infer_haplotype_call(8, 2),
            infer_haplotype_call(1, 1),
        )
    ]
    summary = summarize_analysis([support], calls_by_snv)
    output_stream = StringIO()

    write_analysis_summary(
        summary,
        output_stream,
        eligible_snv_count=1,
        tsv_output="counts.tsv",
        histogram_output="discrepancies.png",
    )

    output = output_stream.getvalue()
    assert "untagged_informative_votes=3" in output
    assert "covered_cells=2" in output
    assert "called_cells=1" in output
    assert "wrote_tsv=counts.tsv" in output
    assert "wrote_histogram=discrepancies.png" in output
