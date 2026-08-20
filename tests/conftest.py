from pathlib import Path

import pysam
import pytest

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def edge_case_bam(tmp_path: Path) -> Path:
    """Convert the reviewable SAM fixture into a temporary indexed BAM."""
    sam_path = FIXTURES / "edge_cases.sam"
    bam_path = tmp_path / "edge_cases.bam"
    with (
        pysam.AlignmentFile(sam_path, "r") as sam_file,
        pysam.AlignmentFile(bam_path, "wb", template=sam_file) as bam_file,
    ):
        for alignment in sam_file:
            bam_file.write(alignment)
    pysam.index(str(bam_path))
    return bam_path
