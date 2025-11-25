"""Module for easy generation of FUSOR AssayedFusion/CategoricalFusion objects"""

import re

from cool_seq_tool.schemas import Assembly, CoordinateType, GenomicTxMetadata, Strand

from fusor.fusor import FUSOR
from fusor.models import AssayedFusion, CategoricalFusion, TranscriptSegmentElement

JUNC_HGVS_PATTERN = re.compile(r"^(?:NM_|NC_)\d+\.\d+:(?:c|g)\.\d+$")


class TranscriptJunctionBuilder:
    """Class for quick FUSOR object generation"""

    def __init__(
        self,
        fusor: FUSOR,
        five_prime_junction: str,
        three_prime_junction: str,
        five_prime_gene: str,
        three_prime_gene: str,
        assayed_fusion: bool = True,
    ) -> None:
        """Initialize TranscriptJunctionBuilder class

        :param fusor: A FUSOR object
        :param five_prime_junction: An HGVS string representation of the
            five prime junction location. This must use a c. or g. coordinate
        :param three_prime_junction: An HGVS string representation of the
            three prime junction location. This must use a c. or g. coordinate
        :param five_prime_gene: The 5' gene partner
        :param three_prime_gene: The 3' gene partner
        :param assayed_fusion: If an `AssayedFusion` object should be created.
            By default, this is set to True.
        :raises ValueError: If ``five_prime_junction`` or
            ``three_prime_junction`` are not described using c. coordinates
        """
        self.fusor = fusor
        for junc in [five_prime_junction, three_prime_junction]:
            if not JUNC_HGVS_PATTERN.match(junc):
                msg = "The fusion junction locations must be described using c. or g. coordinates"
                raise ValueError(msg)
        self.five_prime_junction = five_prime_junction
        self.three_prime_junction = three_prime_junction
        self.five_prime_gene = five_prime_gene
        self.three_prime_gene = three_prime_gene
        self.assayed_fusion = assayed_fusion

    async def _get_cds_start(self, tx: str) -> int:
        """Get CDS start position for transcript

        :param tx: A transcript accession
        :return: The CDS start site for the transcript
        """
        cds_range = await self.fusor.cool_seq_tool.uta_db.get_cds_start_end(tx_ac=tx)
        return int(cds_range[0])

    async def _get_junc_location(
        self, tx: str, pos: int, cds_start: int = 0
    ) -> GenomicTxMetadata:
        """Get junction data from transcript

        :param tx: A transcript accession
        :param pos: The provided junction location
        :param cds_start: The CDS start site
        :return: A `GenomicTxMetadata` object. The returned coordinates will
            always be described on GRCh38
        """
        return await self.fusor.cool_seq_tool.uta_db.get_genomic_tx_data(
            tx_ac=tx,
            pos=((cds_start + pos - 1, cds_start + pos)),
            target_genome_assembly=Assembly.GRCH38,
        )

    async def _extract_junc(self, junc: GenomicTxMetadata) -> int:
        """Extract junction location from GenomicTxMetadata object

        :param junc: A `GenomicTxMetadata` object
        :return: The residue coordinate where the fusion junction occurs
        """
        pos_range = junc.alt_pos_change_range
        return pos_range[1] if junc.strand == Strand.POSITIVE else pos_range[0]

    async def _create_tx_segment(
        self, junc: str, five_prime: bool = True
    ) -> TranscriptSegmentElement:
        """Create TranscriptSegmentElement from junction string

        :param junc: An HGVS string describing the junction location
        :param five_prime: If the 5' prime segment is being created. This is
            set by default to ``True``.
        :return: A `TranscriptSegmentElement` object
        """
        ref_seq = junc.split(":")[0]
        pos = int(junc.split(":")[1].split(".")[1])

        if "c." in junc:
            cds = await self._get_cds_start(ref_seq)
            junc_data = await self._get_junc_location(
                tx=ref_seq, pos=pos, cds_start=cds
            )
            pos = await self._extract_junc(junc=junc_data)
            seg = await self.fusor.transcript_segment_element(
                tx_to_genomic_coords=False,
                transcript=ref_seq,
                genomic_ac=junc_data.alt_ac,
                seg_start_genomic=pos if not five_prime else None,
                seg_end_genomic=pos if five_prime else None,
                coordinate_type=CoordinateType.RESIDUE,
            )
            return seg[0]
        seg = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=ref_seq,
            gene=self.five_prime_gene if five_prime else self.three_prime_gene,
            seg_start_genomic=pos if not five_prime else None,
            seg_end_genomic=pos if five_prime else None,
            coordinate_type=CoordinateType.RESIDUE,
        )
        return seg[0]

    async def build_fusion(
        self,
    ) -> AssayedFusion | CategoricalFusion:
        """Create `AssayedFusion` or `CategoricalFusion` in user-accessible
        way

        :return: An `AssayedFusion` or `CategoricalFusion` object
        """
        five_prime_seg = await self._create_tx_segment(junc=self.five_prime_junction)
        three_prime_seg = await self._create_tx_segment(
            junc=self.three_prime_junction, five_prime=False
        )
        return (
            self.fusor.assayed_fusion(structure=[five_prime_seg, three_prime_seg])
            if self.assayed_fusion
            else self.fusor.categorical_fusion(
                structure=[five_prime_seg, three_prime_seg]
            )
        )
