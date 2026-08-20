"""Core data structures for SNVs and their allele-support counts."""

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class BiallelicSnv:
    """A canonical biallelic SNV with a one-based genomic position."""

    chromosome: str
    position: int
    reference_allele: str
    alternate_allele: str


@dataclass(slots=True)
class HaplotypeSupportCounts:
    """Informative REF and ALT observations for one SNV."""

    snv: BiallelicSnv
    haplotype_1_ref: int = 0
    haplotype_1_alt: int = 0
    haplotype_2_ref: int = 0
    haplotype_2_alt: int = 0
    untagged_ref: int = 0
    untagged_alt: int = 0

    def counts_for_haplotype(self, haplotype: int) -> tuple[int, int]:
        """Return ``(REF, ALT)`` counts for haplotype 1 or 2."""
        if haplotype == 1:
            return self.haplotype_1_ref, self.haplotype_1_alt
        if haplotype == 2:
            return self.haplotype_2_ref, self.haplotype_2_alt
        raise ValueError(f"haplotype must be 1 or 2, not {haplotype}")

    def add_haplotype_vote(self, haplotype: int, allele: str) -> None:
        """Add one REF or ALT vote to haplotype 1 or 2."""
        field_name = f"haplotype_{haplotype}_{allele.lower()}"
        valid_fields = {
            "haplotype_1_ref",
            "haplotype_1_alt",
            "haplotype_2_ref",
            "haplotype_2_alt",
        }
        if field_name not in valid_fields:
            raise ValueError(
                f"unsupported vote: haplotype={haplotype}, allele={allele}"
            )
        setattr(self, field_name, getattr(self, field_name) + 1)

    def add_untagged_vote(self, allele: str) -> None:
        """Record an informative vote from a missing or invalid HP tag."""
        if allele == "REF":
            self.untagged_ref += 1
        elif allele == "ALT":
            self.untagged_alt += 1
        else:
            raise ValueError(f"unsupported allele: {allele}")
