"""Module for easy generation of FUSOR AssayedFusion/CategoricalFusion objects"""

from cool_seq_tool.schemas import (
    AnnotationLayer,
    Assembly,
    CoordinateType,
    GenomicTxMetadata,
    Strand,
)
from ga4gh.core.models import MappableConcept

from fusor.fusor import FUSOR
from fusor.models import (
    AssayedFusion,
    CategoricalFusion,
    GeneElement,
    TranscriptSegmentElement,
)


class TranscriptJunctionBuilder:
    """Class for quick FUSOR object generation"""

    def __init__(
        self,
        fusor: FUSOR,
        five_prime_gene: str,
        three_prime_gene: str,
        five_prime_reference_sequence: str | None = None,
        three_prime_reference_sequence: str | None = None,
        annotation_type: AnnotationLayer | None = None,
        five_prime_junction: int | None = None,
        three_prime_junction: int | None = None,
        assayed_fusion: bool = True,
    ) -> None:
        """Initialize TranscriptJunctionBuilder class

        :param fusor: A FUSOR object
        :param five_prime_gene: The 5' gene partner
        :param three_prime_gene: The 3' gene partner
        :param five_prime_reference_sequence: The 5' prime reference sequence,
            either transcript or genomic. By default, this is set to None
        :param three_prime_reference_sequence: The 5' prime reference sequence,
            either transcript or genomic. By default, this is set to None
        :param annotation_type: The annotation type describing the 5' and 3'
            reference sequences. By default, this is set to None
        :param five_prime_junction: The 5' junction location, described using
            a residue coordinate. By default, this is set to None.
        :param three_prime_junction: The 3' junction location, described using
            a residue coordinate. By default, this is set to None.
        :param assayed_fusion: If an `AssayedFusion` object should be created.
            By default, this is set to True.
        :raises ValueError: If ``five_prime_junction`` or
            ``three_prime_junction`` are not described using c. coordinates
        """
        self.fusor = fusor

        # Validate only c. or g. annotation type is provided
        if annotation_type == AnnotationLayer.PROTEIN:
            msg = "Only c. or g. RefSeq accessions are supported"
            raise ValueError(msg)

        # Validate that five_prime_junction and three_prime_junction are
        # provided if their accessions are provided
        if five_prime_reference_sequence and not five_prime_junction:
            msg = "Please provide 5' junction location"
            raise ValueError(msg)
        if three_prime_reference_sequence and not three_prime_junction:
            msg = "Please provide 3' junction location"
            raise ValueError(msg)

        self.five_prime_gene = five_prime_gene
        self.three_prime_gene = three_prime_gene
        self.five_prime_reference_sequence = five_prime_reference_sequence
        self.three_prime_reference_sequence = three_prime_reference_sequence
        self.annotation_type = annotation_type
        self.five_prime_junction = five_prime_junction
        self.three_prime_junction = three_prime_junction
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
        self, ref_seq: str, pos: int, five_prime: bool = True
    ) -> TranscriptSegmentElement:
        """Create TranscriptSegmentElement from junction string

        :param ref_seq: The reference sequence for the fusion partner
        :param pos: The fusion junction location
        :param five_prime: If the 5' prime segment is being created. This is
            set by default to ``True``.
        :return: A `TranscriptSegmentElement` object
        """
        if self.annotation_type == AnnotationLayer.CDNA:
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

    async def _process_partner(
        self,
        gene: str,
        ref_seq: str | None = None,
        pos: int | None = None,
        five_prime: bool = True,
    ) -> GeneElement | TranscriptSegmentElement:
        """Process fusion partner input

        :param gene: The gene symbol for the fusion partner
        :param ref_seq: The chromosomal or transcript reference sequence.
            By default, this is set to None
        :param pos: The fusion junction location. By default, this is set to
            None
        :param five_prime: If the 5' prime partner is being examined. By
            default, this is set to True
        """
        if not ref_seq:
            gene_obj = self.fusor.gene_element(gene=gene)[0]
            if gene_obj:
                return gene_obj
            if not gene_obj:
                return GeneElement(gene=MappableConcept(name=gene, conceptType="Gene"))
        return await self._create_tx_segment(
            ref_seq=ref_seq, pos=pos, five_prime=five_prime
        )

    async def build_fusion(
        self,
    ) -> AssayedFusion | CategoricalFusion:
        """Create `AssayedFusion` or `CategoricalFusion` in user-accessible
        way

        :return: An `AssayedFusion` or `CategoricalFusion` object
        """
        five_prime_seg = await self._process_partner(
            gene=self.five_prime_gene,
            ref_seq=self.five_prime_reference_sequence,
            pos=self.five_prime_junction,
            five_prime=True,
        )
        three_prime_seg = await self._process_partner(
            gene=self.three_prime_gene,
            ref_seq=self.three_prime_reference_sequence,
            pos=self.three_prime_junction,
            five_prime=False,
        )
        return (
            self.fusor.assayed_fusion(structure=[five_prime_seg, three_prime_seg])
            if self.assayed_fusion
            else self.fusor.categorical_fusion(
                structure=[five_prime_seg, three_prime_seg]
            )
        )
