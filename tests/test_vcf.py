from pathlib import Path

from hap_counter.vcf import load_biallelic_snvs

TEST_DATA = Path(__file__).parents[1] / "test_data"
PLAIN_VCF = TEST_DATA / "giab_2023.05.hg002.wf_snp.chr16_28000000_29000000.vcf"
COMPRESSED_VCF = PLAIN_VCF.with_suffix(".vcf.gz")


def test_plain_and_compressed_vcfs_produce_the_same_snvs() -> None:
    plain = load_biallelic_snvs(PLAIN_VCF)
    compressed = load_biallelic_snvs(COMPRESSED_VCF)

    assert plain == compressed
    assert len(plain) == 939
    assert plain[0].chromosome == "chr16"
    assert plain[0].position == 28_001_381


def test_fixture_filters_noncanonical_and_non_snv_records() -> None:
    fixture_vcf = Path(__file__).parent / "fixtures" / "variants.vcf"

    snvs = load_biallelic_snvs(fixture_vcf)

    assert [snv.position for snv in snvs] == [101, 201, 301, 401, 501]
