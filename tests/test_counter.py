from io import StringIO
from pathlib import Path

import pytest

from hap_counter.counter import ReadFilterConfig, count_haplotype_support
from hap_counter.models import BiallelicSnv
from hap_counter.output import write_support_tsv
from hap_counter.vcf import load_biallelic_snvs

FIXTURES = Path(__file__).parent / "fixtures"
TEST_DATA = Path(__file__).parents[1] / "test_data"


def test_counts_match_reviewable_expected_output(edge_case_bam: Path) -> None:
    snvs = load_biallelic_snvs(FIXTURES / "variants.vcf")
    support_counts = count_haplotype_support(edge_case_bam, snvs)
    output_stream = StringIO()

    write_support_tsv(support_counts, output_stream)

    assert output_stream.getvalue() == (
        FIXTURES / "expected_counts.tsv"
    ).read_text()


def test_cigar_edge_cases_are_mapped_to_the_correct_base(
    edge_case_bam: Path,
) -> None:
    snvs = load_biallelic_snvs(FIXTURES / "variants.vcf")
    support_by_position = {
        item.snv.position: item
        for item in count_haplotype_support(edge_case_bam, snvs)
    }

    assert support_by_position[201].counts_for_haplotype(1) == (0, 1)
    assert support_by_position[301].counts_for_haplotype(1) == (0, 0)
    assert support_by_position[401].counts_for_haplotype(1) == (0, 1)


def test_quality_thresholds_are_configurable(edge_case_bam: Path) -> None:
    support = count_haplotype_support(
        edge_case_bam,
        [BiallelicSnv("chr1", 101, "C", "G")],
        ReadFilterConfig(minimum_mapping_quality=10, minimum_base_quality=10),
    )[0]

    assert (support.haplotype_1_ref, support.haplotype_1_alt) == (1, 1)


def test_untagged_and_invalid_hp_votes_are_diagnostic(edge_case_bam: Path) -> None:
    support = count_haplotype_support(
        edge_case_bam, [BiallelicSnv("chr1", 101, "C", "G")]
    )[0]

    assert (support.untagged_ref, support.untagged_alt) == (0, 2)


def test_missing_bam_contig_is_reported(edge_case_bam: Path) -> None:
    with pytest.raises(ValueError, match="chr2"):
        count_haplotype_support(
            edge_case_bam, [BiallelicSnv("chr2", 101, "C", "G")]
        )


def test_supplied_bam_counts_first_fixture_snv() -> None:
    bam_path = TEST_DATA / (
        "giab_2023.05.hg002.haplotagged."
        "chr16_28000000_29000000.processed.30x.bam"
    )

    support = count_haplotype_support(
        bam_path, [BiallelicSnv("chr16", 28_001_381, "G", "A")]
    )[0]

    assert (support.haplotype_1_ref, support.haplotype_1_alt) == (11, 4)
    assert (support.haplotype_2_ref, support.haplotype_2_alt) == (3, 12)
