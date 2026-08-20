"""VCF input handling."""

from pathlib import Path

import pysam

from hap_counter.models import BiallelicSnv

CANONICAL_BASES = frozenset("ACGT")


def load_biallelic_snvs(vcf_path: str | Path) -> list[BiallelicSnv]:
    """Load canonical biallelic SNVs from a VCF in file order."""
    snvs: list[BiallelicSnv] = []
    with pysam.VariantFile(str(vcf_path)) as variant_file:
        for variant_record in variant_file:
            alternate_alleles = variant_record.alts or ()
            if len(alternate_alleles) != 1:
                continue

            reference_allele = variant_record.ref.upper()
            alternate_allele = alternate_alleles[0].upper()
            if (
                reference_allele not in CANONICAL_BASES
                or alternate_allele not in CANONICAL_BASES
                or reference_allele == alternate_allele
            ):
                continue

            snvs.append(
                BiallelicSnv(
                    chromosome=variant_record.contig,
                    position=variant_record.pos,
                    reference_allele=reference_allele,
                    alternate_allele=alternate_allele,
                )
            )
    return snvs
