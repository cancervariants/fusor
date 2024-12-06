"""Module for translating output from fusion detection methods to FUSOR AssayedFusion
objects
"""

import logging
from enum import Enum

import polars as pl
from cool_seq_tool.schemas import Assembly

from fusor.fusor import FUSOR
from fusor.models import (
    Assay,
    AssayedFusion,
    CausativeEvent,
    EventType,
    GeneElement,
    TranscriptSegmentElement,
)

_logger = logging.getLogger(__name__)


class Caller(str, Enum):
    """Define different supported callers"""

    JAFFA = "JAFFA"
    STAR_FUSION = "STAR-Fusion"
    FUSION_CATCHER = "FusionCatcher"
    FUSION_MAP = "FusionMap"
    ARRIBA = "Arriba"
    CICERO = "CICERO"
    MAPSPLICE = "MapSplice"
    ENFUSION = "EnFusion"
    GENIE = "GENIE"


class Translator:
    """Class for translating outputs from different fusion detection algorithms
    to FUSOR AssayedFusion objects
    """

    def __init__(self, fusor: FUSOR) -> None:
        """Initialize Translator class

        :param fusor: A FUSOR instance
        """
        self.fusor = fusor

    def _format_fusion(
        self,
        gene_5prime: GeneElement,
        gene_3prime: GeneElement,
        tr_5prime: TranscriptSegmentElement | None = None,
        tr_3prime: TranscriptSegmentElement | None = None,
        ce: CausativeEvent | None = None,
        rf: bool | None = None,
        assay: Assay | None = None,
    ) -> AssayedFusion:
        """Format classes to create AssayedFusion objects

        :param gene_5prime: 5'prime GeneElement
        :param gene_3prime: 3'prime GeneElement
        :param tr_5prime: 5'prime TranscriptSegmentElement
        :param tr_3prime: 3'prime TranscriptSegmentElement
        :param ce: CausativeEvent
        :param rf: A boolean indicating if the reading frame is preserved
        :param assay: Assay
        :return AssayedFusion object
        """
        params = {}
        if not tr_5prime[0] and not tr_3prime[0]:
            params["structure"] = [gene_5prime, gene_3prime]
        elif tr_5prime[0] and not tr_3prime[0]:
            params["structure"] = [tr_5prime[0], gene_3prime]
        elif not tr_5prime[0] and tr_3prime[0]:
            params["structure"] = [gene_5prime, tr_3prime[0]]
        else:
            params["structure"] = [tr_5prime[0], tr_3prime[0]]

        if ce:
            params["causativeEvent"] = ce
        if rf:
            params["readingFramePreserved"] = rf
        if assay:
            params["assay"] = assay
        return AssayedFusion(**params)

    def _get_causative_event(
        self, chrom1: str, chrom2: str, descr: str | None = None
    ) -> CausativeEvent | None:
        """Infer Causative Event. Currently restricted to rearrangements

        :param chrom1: The chromosome for the 5' partner
        :param chrom2: The chromosome for the 3' partner
        :param descr: An annotation describing the fusion event. This input is supplied to the eventDescription CausativeEvent attribute.
        :return: A CausativeEvent object if construction is successful
        """
        if descr and "rearrangement" in descr:
            return CausativeEvent(
                eventType=EventType("rearrangement"), eventDescription=descr
            )
        if chrom1 != chrom2:
            return CausativeEvent(eventType=EventType("rearrangement"))
        return None

    def _get_gene_element_unnormalized(self, symbol: str) -> GeneElement:
        """Return GeneElement when gene symbol cannot be normalized

        :param symbol: A gene symbol for a fusion partner
        :return: A GeneElement object
        """
        return GeneElement(
            gene={
                "id": f"gene:{symbol}",
                "label": symbol,
                "type": "Gene",
            },
        )

    def _get_gene_element(self, genes: str, caller: Caller) -> GeneElement:
        """Return a GeneElement given an individual/list of gene symbols and a
        fusion detection algorithm

        :param genes: A gene symbol or list of gene symbols, separated by columns
        :param caller: The examined fusion detection algorithm
        :return A GeneElement object
        """
        if "," not in genes or caller != caller.ARRIBA:
            ge = self.fusor.gene_element(gene=genes)
            return ge[0] if ge[0] else self._get_gene_element_unnormalized(genes)

        genes = genes.split(",")
        dists = []
        for gene in genes:
            start, end = gene.rfind("("), gene.rfind(")")
            dists.append(int(gene[start + 1 : end]))
        gene = (
            genes[0].split("(")[0] if dists[0] <= dists[1] else genes[1].split("(")[0]
        )
        ge = self.fusor.gene_element(gene=gene)
        return ge[0] if ge[0] else self._get_gene_element_unnormalized(gene)

    def _get_genomic_ac(self, chrom: str, build: Assembly) -> str:
        """Return a RefSeq genomic accession given a chromosome and a reference build

        :param chrom: A chromosome number
        :param build: The assembly, either GRCh37 or GRCh38
        :return: The corresponding refseq genomic accession
        :raise ValueError: if unable to retrieve genomic accession
        """
        sr = self.fusor.cool_seq_tool.seqrepo_access
        alias_list, errors = sr.translate_identifier(
            f"{build}:{chrom}", target_namespaces="refseq"
        )
        if errors:
            statement = f"Genomic accession for {chrom} could not be retrieved"
            _logger.error(statement)
            raise ValueError
        return alias_list[0].split(":")[1]

    async def from_jaffa(
        self,
        fusion_genes: str,
        chrom1: str,
        base1: int,
        chrom2: str,
        base2: int,
        rearrangement: bool,
        classification: str,
        inframe: bool,
        rb: Assembly,
    ) -> AssayedFusion | None:
        """Parse JAFFA fusion output to create AssayedFusion object

        :param fusion_genes: A string containing the two fusion partners
        :param chrom1: The chromosome indicated in the chrom1 column
        :param base1: The genomic position indicated in the base1 column
        :param chrom2: The chromosome indicated in the chrom2 column
        :param base2: The genomic position indicated in the base2 column
        :param rearrangement: A boolean indicating if a rearrangement occured
        :param classification: The classification associated with the called fusion
        :param inframe: A boolean indicating if the fusion occurred in-frame
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        genes = fusion_genes.split(":")
        gene_5prime_element = self._get_gene_element(genes[0], Caller.JAFFA)
        gene_3prime_element = self._get_gene_element(genes[1], Caller.JAFFA)
        gene_5prime = gene_5prime_element.gene.label
        gene_3prime = gene_3prime_element.gene.label

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(chrom1, rb),
            seg_end_genomic=base1,
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(chrom2, rb),
            seg_start_genomic=base2,
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        if rearrangement:
            ce = CausativeEvent(
                eventType=EventType("rearrangement"),
                eventDescription=classification,
            )
        else:
            ce = None

        return self._format_fusion(
            gene_5prime_element, gene_3prime_element, tr_5prime, tr_3prime, ce, inframe
        )

    async def from_star_fusion(
        self,
        left_gene: str,
        right_gene: str,
        left_breakpoint: str,
        right_breakpoint: str,
        annots: str,
        rb: Assembly,
    ) -> AssayedFusion:
        """Parse STAR-Fusion output to create AssayedFusion object

        :param left_gene: The gene indicated in the LeftGene column
        :param right_gene: The gene indicated in the RightGene column
        :param left_breakpoint: The gene indicated in the LeftBreakpoint column
        :param right_breakpoint: The gene indicated in the RightBreakpoint column
        :param annots: The annotations associated with the fusion
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = left_gene.split("^")[0]
        gene2 = right_gene.split("^")[0]
        gene_5prime_element = self._get_gene_element(gene1, Caller.STAR_FUSION)
        gene_3prime_element = self._get_gene_element(gene2, Caller.STAR_FUSION)
        gene_5prime = gene_5prime_element.gene.label
        gene_3prime = gene_3prime_element.gene.label

        five_prime = left_breakpoint.split(":")
        three_prime = right_breakpoint.split(":")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(five_prime[0], rb),
            seg_end_genomic=int(five_prime[1]),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(three_prime[0], rb),
            seg_start_genomic=int(three_prime[1]),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(five_prime[0], three_prime[0], ",".join(annots))
        return self._format_fusion(
            gene_5prime_element, gene_3prime_element, tr_5prime, tr_3prime, ce
        )

    async def from_fusion_catcher(
        self,
        five_prime_partner: str,
        three_prime_partner: str,
        five_prime_fusion_point: str,
        three_prime_fusion_point: str,
        predicted_effect: str,
        rb: Assembly,
    ) -> AssayedFusion:
        """Parse FusionCatcher output to create AssayedFusion object

        :param five_prime_partner: Gene symbol for the 5' fusion partner
        :param three_prime_partner: Gene symbol for the 3' fusion partner
        :param five_prime_fusion_point: Chromosomal position for the 5' end of the
        fusion junction. This coordinate is 1-based
        :param three_prime_fusion_point:  Chromosomal position for the 3' end of the
        fusion junction. This coordinate is 1-based
        :param predicted_effect: The predicted effect of the fusion event, created
        using annotation from the Ensembl database
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene_5prime_element = self._get_gene_element(
            five_prime_partner, Caller.FUSION_CATCHER
        )
        gene_3prime_element = self._get_gene_element(
            three_prime_partner, Caller.FUSION_CATCHER
        )

        five_prime = five_prime_fusion_point.split(":")
        three_prime = three_prime_fusion_point.split(":")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(five_prime[0], rb),
            seg_end_genomic=int(five_prime[1]),
            gene=gene_5prime_element.gene.label,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(three_prime[0], rb),
            seg_start_genomic=int(three_prime[1]),
            gene=gene_3prime_element.gene.label,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(five_prime[0], three_prime[0], predicted_effect)
        return self._format_fusion(
            gene_5prime_element, gene_3prime_element, tr_5prime, tr_3prime, ce
        )

    async def from_fusion_map(
        self, fmap_row: pl.DataFrame, rb: Assembly
    ) -> AssayedFusion:
        """Parse FusionMap output to create FUSOR AssayedFusion object

        :param fmap_row: A row of FusionMap output
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = fmap_row.get_column("KnownGene1").item()
        gene2 = fmap_row.get_column("KnownGene2").item()
        gene_5prime = self._get_gene_element(gene1, "fusion_map").gene.label
        gene_3prime = self._get_gene_element(gene2, "fusion_map").gene.label

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(
                fmap_row.get_column("Chromosome1").item(), rb
            ),
            seg_end_genomic=int(fmap_row.get_column("Position1").item()),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(
                fmap_row.get_column("Chromosome2").item(), rb
            ),
            seg_start_genomic=int(fmap_row.get_column("Position2").item()),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        # Combine columns to create fusion annotation string"
        descr = (
            fmap_row.get_column("FusionGene").item()
            + ","
            + fmap_row.get_column("SplicePatternClass").item()
            + ","
            + fmap_row.get_column("FrameShiftClass").item()
        )
        ce = self._get_causative_event(
            fmap_row.get_column("Chromosome1").item(),
            fmap_row.get_column("Chromosome2").item(),
            descr,
        )
        rf = bool(fmap_row.get_column("FrameShiftClass").item() == "InFrame")
        return self._format_fusion(
            gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce, rf
        )

    async def from_arriba(
        self,
        gene1: str,
        gene2: str,
        strand1: str,
        strand2: str,
        breakpoint1: str,
        breakpoint2: str,
        event: str,
        confidence: str,
        direction1: str,
        direction2: str,
        rf: str,
        rb: Assembly,
    ) -> AssayedFusion:
        """Parse Arriba output to create AssayedFusion object

        :param gene1: The 5' gene fusion partner
        :param gene2: The 3' gene fusion partner
        :param strand1: The strand information for the 5' gene fusion partner
        :param strand2: The strand information for the 3' gene fusion partner
        :param breakpoint1: The chromosome and breakpoint for gene1
        :param breakpoint2: The chromosome and breakpoint for gene2
        :param event: An inference about the type of fusion event
        :param confidence: A metric describing the confidence of the fusion prediction
        :param direction1: A description that indicates if the transcript segment
            starts or ends at breakpoint1
        :param direction2: A description that indicates if the transcript segment
            starts or ends at breakpoint2
        :param rf: A description if the reading frame is preserved for the fusion
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        # Arriba reports two gene symbols if a breakpoint occurs in an intergenic
        # space. We select the gene symbol with the smallest distance from the
        # breakpoint.
        gene_5prime_element = self._get_gene_element(gene1, "arriba")
        gene_3prime_element = self._get_gene_element(gene2, "arriba")
        gene_5prime = gene_5prime_element.gene.label
        gene_3prime = gene_3prime_element.gene.label

        strand1 = strand1.split("/")[1]  # Determine strand that is transcribed
        strand2 = strand2.split("/")[1]  # Determine strand that is transcribed
        if strand1 == "+":
            gene1_seg_start = direction1 == "upstream"
        else:
            gene1_seg_start = direction1 == "downstream"
        if strand2 == "+":
            gene2_seg_start = direction2 == "upstream"
        else:
            gene2_seg_start = direction2 == "downstream"

        breakpoint1 = breakpoint1.split(":")
        breakpoint2 = breakpoint2.split(":")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(breakpoint1[0], rb),
            seg_start_genomic=int(breakpoint1[1]) if gene1_seg_start else None,
            seg_end_genomic=int(breakpoint1[1]) if not gene1_seg_start else None,
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(breakpoint2[0], rb),
            seg_start_genomic=int(breakpoint2[1]) if gene2_seg_start else None,
            seg_end_genomic=int(breakpoint2[1]) if not gene2_seg_start else None,
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        ce = (
            CausativeEvent(
                eventType=EventType("read-through"),
                eventDescription=confidence,
            )
            if "read_through" in event
            else CausativeEvent(
                eventType=EventType("rearrangement"),
                eventDescription=confidence,
            )
        )
        rf = bool(rf == "in-frame") if rf != "." else None
        return self._format_fusion(
            gene_5prime_element, gene_3prime_element, tr_5prime, tr_3prime, ce, rf
        )

    async def from_cicero(
        self,
        gene_5prime: str,
        gene_3prime: str,
        chr_5prime: str,
        chr_3prime: str,
        pos_5prime: int,
        pos_3prime: int,
        sv_ort: str,
        event_type: str,
        rb: Assembly,
    ) -> AssayedFusion | str:
        """Parse CICERO output to create AssayedFusion object

        :param gene_5prime: The gene symbol for the 5' partner
        :param gene_3prime: The gene symbol for the 3' partner
        :param chr_5prime: The chromosome for the 5' partner
        :param chr_3prime: The chromosome for the 3' partner
        :param pos_5prime: The genomic breakpoint for the 5' partner
        :param pos_3prime: The genomic breakpoint for the 3' partner
        :param sv_ort: Whether the mapping orientation of assembled contig (driven by
            structural variation) has confident biological meaning
        :param event_type: The structural variation event that created the called fusion
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        # Check if gene symbols have valid formatting. CICERO can output two or more
        # gene symbols for `gene_5prime` or `gene_3prime`, which are separated by a comma. As
        # there is not a precise way to resolve this ambiguity, we do not process
        # these events
        if "," in gene_5prime or "," in gene_3prime:
            msg = "Ambiguous gene symbols are reported by CICERO for at least one of the fusion partners"
            _logger.warning(msg)
            return msg

        # Check CICERO annotation regarding the confidence that the called fusion
        # has biological meaning
        if sv_ort != ">":
            msg = "CICERO annotation indicates that this event does not have confident biological meaning"
            _logger.warning(msg)
            return msg

        gene_5prime_element = self._get_gene_element(gene_5prime, "cicero")
        gene_3prime_element = self._get_gene_element(gene_3prime, "cicero")
        gene_5prime = gene_5prime_element.gene.label
        gene_3prime = gene_3prime_element.gene.label

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(chr_5prime, rb),
            seg_end_genomic=pos_5prime,
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(chr_3prime, rb),
            seg_start_genomic=pos_3prime,
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        if event_type == "read_through":
            ce = CausativeEvent(
                eventType=EventType("read-through"),
                eventDescription=event_type,
            )
        else:
            ce = CausativeEvent(
                eventType=EventType("rearrangement"),
                eventDescription=event_type,
            )
        return self._format_fusion(
            gene_5prime_element,
            gene_3prime_element,
            tr_5prime,
            tr_3prime,
            ce,
        )

    async def from_mapsplice(
        self, mapsplice_row: pl.DataFrame, rb: Assembly
    ) -> AssayedFusion:
        """Parse MapSplice output to create AssayedFusion object

        :param mapsplice_row: A row of MapSplice output
        :param rb: The reference build used to call the fusion
        :retun: An AssayedFusion object, if construction is successful
        """
        gene1 = mapsplice_row[60].strip(",")
        gene2 = mapsplice_row[61].strip(",")
        gene_5prime_element = self._get_gene_element(gene1, "mapsplice")
        gene_3prime_element = self._get_gene_element(gene2, "mapsplice")
        gene_5prime = gene_5prime_element.gene.label
        gene_3prime = gene_3prime_element.gene.label

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(mapsplice_row[0].split("~")[0], rb),
            seg_end_genomic=int(mapsplice_row[1]),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(mapsplice_row[0].split("~")[1], rb),
            seg_start_genomic=int(mapsplice_row[2]),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(
            mapsplice_row[0].split("~")[0], mapsplice_row[0].split("~")[1]
        )
        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)

    async def from_enfusion(
        self,
        gene_5prime: str,
        gene_3prime: str,
        chr_5prime: int,
        chr_3prime: int,
        break_5prime: int,
        break_3prime: int,
        rb: Assembly,
    ) -> AssayedFusion:
        """Parse EnFusion output to create AssayedFusion object

        :param gene_5prime: The 5' gene fusion partner
        :param gene_3prime: The 3' gene fusion partner
        :param chr_5prime: The 5' gene fusion partner chromosome
        :param chr_3prime: The 3' gene fusion partner chromosome
        :param break_5prime: The 5' gene fusion partner genomic breakpoint
        :param break_3prime: The 3' gene fusion partner genomic breakpoint
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene_5prime_element = self._get_gene_element(gene_5prime, "enfusion")
        gene_3prime_element = self._get_gene_element(gene_3prime, "enfusion")
        gene_5prime = gene_5prime_element.gene.label
        gene_3prime = gene_3prime_element.gene.label

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(chr_5prime, rb),
            seg_end_genomic=break_5prime,
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(chr_3prime, rb),
            seg_start_genomic=break_3prime,
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(
            chr_5prime,
            chr_3prime,
        )
        return self._format_fusion(
            gene_5prime_element, gene_3prime_element, tr_5prime, tr_3prime, ce
        )

    async def from_genie(
        self,
        site1_hugo: str,
        site2_hugo: str,
        site1_chrom: int,
        site2_chrom: int,
        site1_pos: int,
        site2_pos: int,
        annot: str,
        reading_frame: str,
        rb: Assembly,
    ) -> AssayedFusion:
        """Parse GENIE output to create AssayedFusion object

        :param site1_hugo: The HUGO symbol reported at site 1
        :param site2_hugo: The HUGO symbol reported at site 2
        :param site1_chrom: The chromosome reported at site 1
        :param site2_chrom: The chromosome reported at site 2
        :param site1_pos: The breakpoint reported at site 1
        :param site2_pos: The breakpoint reported at site 2
        :param annot: The annotation for the fusion event
        :param reading_frame: The reading frame status of the fusion
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene_5prime_element = self._get_gene_element(site1_hugo, "genie")
        gene_3prime_element = self._get_gene_element(site2_hugo, "genie")
        gene_5prime = gene_5prime_element.gene.label
        gene_3prime = gene_3prime_element.gene.label

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(site1_chrom, rb),
            seg_end_genomic=site1_pos,
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(site2_chrom, rb),
            seg_start_genomic=site2_pos,
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(
            site1_chrom,
            site2_chrom,
            annot,
        )
        rf = bool(reading_frame == "in frame")
        return self._format_fusion(
            gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce, rf
        )
