"""Module for translating output from fusion detection methods to FUSOR AssayedFusion
objects
"""

import logging

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
        :param Optional[tr_5prime]: 5'prime TranscriptSegmentElement
        :param Optional[tr_3prime]: 3'prime TranscriptSegmentElement
        :param Optional[ce]: CausativeEvent
        :param Optional[rf]: A boolean indicating if the reading frame is preserved
        :param Optional[assay]: Assay
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
        return AssayedFusion(**params), None

    def _get_causative_event(
        self, chrom1: str, chrom2: str, descr: str | None = None
    ) -> CausativeEvent:
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
        return (
            GeneElement(
                type="GeneElement",
                gene={
                    "id": f"gene:{symbol}",
                    "label": symbol,
                    "type": "Gene",
                },
            ),
            None,
        )

    def _get_gene_element(self, genelist: str, caller: str) -> GeneElement:
        """Return a GeneElement given an individual/list of gene symbols and a
        fusion detection algorithm

        :param genelist: A gene symbol or list of gene symbols
        :param caller: The examined fusion detection algorithm
        :return A GeneElement object
        """
        if "," not in genelist or caller != "arriba":
            ge = self.fusor.gene_element(gene=genelist)
            return ge if ge[0] else self._get_gene_element_unnormalized(genelist)

        genes = genelist.split(",")
        dists = []
        for gene in genes:
            start, end = gene.rfind("("), gene.rfind(")")
            dists.append(int(gene[start + 1 : end]))
        gene = (
            genes[0].split("(")[0] if dists[0] <= dists[1] else genes[1].split("(")[0]
        )
        ge = self.fusor.gene_element(gene=gene)
        return ge if ge[0] else self._get_gene_element_unnormalized(gene)

    def _get_genomic_ac(self, chrom: str, build: Assembly) -> str:
        """Return a RefSeq genomic accession given a chromosome and a reference build
        :param chrom: A chromosome number

        :param build: The assembly, either GRCh37 or GRCh38
        :return The corresponding refseq genomic accession
        """
        sr = self.fusor.cool_seq_tool.seqrepo_access
        if build == Assembly.GRCH37:
            alias_list, errors = sr.translate_identifier(f"GRCh37:{chrom}")
        else:
            alias_list, errors = sr.translate_identifier(f"GRCh38:{chrom}")
        if errors:
            _logger.error("Genomic accession for {chrom} could not be retrieved")
            raise ValueError
        for alias in alias_list:
            if alias.startswith("refseq:"):
                genomic_ac = alias.split(":")[1]
        return genomic_ac

    async def from_jaffa(
        self, jaffa_row: pl.DataFrame, rb: Assembly
    ) -> AssayedFusion | None:
        """Parse JAFFA fusion output to create AssayedFusion object

        :param jaffa_row: A row of JAFFA output
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        genes = jaffa_row.get_column("fusion genes").item().split(":")
        gene_5prime_element = self._get_gene_element(genes[0], "jaffa")[0]
        gene_3prime_element = self._get_gene_element(genes[1], "jaffa")[0]
        gene_5prime = gene_5prime_element.gene.label
        gene_3prime = gene_3prime_element.gene.label

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(jaffa_row.get_column("chrom1").item(), rb),
            seg_end_genomic=int(jaffa_row.get_column("base1").item()),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(jaffa_row.get_column("chrom2").item(), rb),
            seg_start_genomic=int(jaffa_row.get_column("base2").item()),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        if jaffa_row.get_column("rearrangement").item() == "TRUE":
            ce = CausativeEvent(
                eventType=EventType("rearrangement"),
                eventDescription=jaffa_row.get_column("classification").item(),
            )
        else:
            ce = None

        rf = bool(jaffa_row.get_column("inframe").item() == "TRUE")
        return self._format_fusion(
            gene_5prime_element, gene_3prime_element, tr_5prime, tr_3prime, ce, rf
        )

    async def from_star_fusion(
        self, sf_row: pl.DataFrame, rb: Assembly
    ) -> AssayedFusion:
        """Parse STAR-Fusion output to create AssayedFusion object

        :param sf_row: A row of STAR-Fusion output
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = sf_row.get_column("LeftGene").item().split("^")[0]
        gene2 = sf_row.get_column("RightGene").item().split("^")[0]
        gene_5prime_element = self._get_gene_element(gene1, "star_fusion")[0]
        gene_3prime_element = self._get_gene_element(gene2, "star_fusion")[0]
        gene_5prime = gene_5prime_element.gene.label
        gene_3prime = gene_3prime_element.gene.label

        five_prime = sf_row.get_column("LeftBreakpoint").item().split(":")
        three_prime = sf_row.get_column("RightBreakpoint").item().split(":")

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

        ce = self._get_causative_event(
            five_prime[0], three_prime[0], ",".join(sf_row.get_column("annots").item())
        )
        return self._format_fusion(
            gene_5prime_element, gene_3prime_element, tr_5prime, tr_3prime, ce
        )

    async def from_fusion_catcher(
        self, fc_row: pl.DataFrame, rb: Assembly
    ) -> AssayedFusion:
        """Parse FusionCatcher output to create AssayedFusion object

        :param fc_row: A row of FusionCatcher output
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = fc_row.get_column("Gene_1_symbol(5end_fusion_partner)").item()
        gene2 = fc_row.get_column("Gene_2_symbol(3end_fusion_partner)").item()
        gene_5prime = self._get_gene_element(gene1, "fusion_catcher")[0].gene.label
        gene_3prime = self._get_gene_element(gene2, "fusion_catcher")[0].gene.label

        five_prime = (
            fc_row.get_column("Fusion_point_for_gene_1(5end_fusion_partner)")
            .item()
            .split(":")
        )
        three_prime = (
            fc_row.get_column("Fusion_point_for_gene_2(3end_fusion_partner)")
            .item()
            .split(":")
        )

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

        ce = self._get_causative_event(
            five_prime[0], three_prime[0], fc_row.get_column("Predicted_effect").item()
        )
        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)

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
        gene_5prime = self._get_gene_element(gene1, "fusion_map")[0].gene.label
        gene_3prime = self._get_gene_element(gene2, "fusion_map")[0].gene.label

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
        self, arriba_row: pl.DataFrame, rb: Assembly
    ) -> AssayedFusion:
        """Parse Arriba output to create AssayedFusion object

        :param arriba_row: A row of Arriba output
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = arriba_row.get_column("#gene1").item()
        gene2 = arriba_row.get_column("gene2").item()

        # Arriba reports two gene symbols if a breakpoint occurs in an intergenic
        # space. We select the gene symbol with the smallest distance from the
        # breakpoint.
        gene_5prime = self._get_gene_element(gene1, "arriba")[0].gene.label
        gene_3prime = self._get_gene_element(gene2, "arriba")[0].gene.label

        breakpoint1 = arriba_row.get_column("breakpoint1").item().split(":")
        breakpoint2 = arriba_row.get_column("breakpoint2").item().split(":")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(breakpoint1[0], rb),
            seg_end_genomic=int(breakpoint1[1]),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(breakpoint2[0], rb),
            seg_start_genomic=int(breakpoint2[1]),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        ce = (
            CausativeEvent(
                eventType=EventType("read-through"),
                eventDescription=arriba_row.get_column("confidence").item(),
            )
            if "read_through" in arriba_row["type"]
            else CausativeEvent(
                eventType=EventType("rearrangement"),
                eventDescription=arriba_row.get_column("confidence").item(),
            )
        )
        rf = bool(arriba_row.get_column("reading_frame").item() == "in-frame")
        return self._format_fusion(
            gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce, rf
        )

    async def from_cicero(
        self, cicero_row: pl.DataFrame, rb: Assembly
    ) -> AssayedFusion:
        """Parse CICERO output to create AssayedFusion object

        :param cicero_row: A row of CICERO output
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = cicero_row.get_column("geneA").item()
        gene2 = cicero_row.get_column("geneB").item()
        gene_5prime = self._get_gene_element(gene1, "cicero")[0].gene.label
        gene_3prime = self._get_gene_element(gene2, "cicero")[0].gene.label

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(cicero_row.get_column("chrA").item(), rb),
            seg_end_genomic=int(cicero_row.get_column("posA").item()),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(cicero_row.get_column("chrB").item(), rb),
            seg_start_genomic=int(cicero_row.get_column("posB").item()),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        if cicero_row.get_column("type").item() == "read_through":
            ce = CausativeEvent(
                eventType=EventType("read-through"),
                eventDescription=cicero_row.get_column("type").item(),
            )
        else:
            ce = CausativeEvent(
                eventType=EventType("rearrangement"),
                eventDescription=cicero_row.get_column("type").item(),
            )
        return self._format_fusion(
            gene_5prime,
            gene_3prime,
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
        gene_5prime = self._get_gene_element(gene1, "mapsplice")[0].gene.label
        gene_3prime = self._get_gene_element(gene2, "mapsplice")[0].gene.label

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
        self, enfusion_row: pl.DataFrame, rb: Assembly
    ) -> AssayedFusion:
        """Parse EnFusion output to create AssayedFusion object

        :param enfusion_row: A row of EnFusion output
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = enfusion_row.get_column("Gene1").item()
        gene2 = enfusion_row.get_column("Gene2").item()
        gene_5prime = self._get_gene_element(gene1, "enfusion")[0].gene.label
        gene_3prime = self._get_gene_element(gene2, "enfusion")[0].gene.label

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(enfusion_row.get_column("Chr1").item(), rb),
            seg_end_genomic=int(enfusion_row.get_column("Break1").item()),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(enfusion_row.get_column("Chr2").item(), rb),
            seg_start_genomic=int(enfusion_row.get_column("Break2").item()),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(
            enfusion_row.get_column("Chr1").item(),
            enfusion_row.get_column("Chr2").item(),
        )
        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)

    async def from_genie(self, genie_row: pl.DataFrame, rb: Assembly) -> AssayedFusion:
        """Parse GENIE output to create AssayedFusion object

        :param genie_row: A row of EnFusion output
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = genie_row.get_column("Site1_Hugo_Symbol").item()
        gene2 = genie_row.get_column("Site2_Hugo_Symbol").item()
        gene_5prime = self._get_gene_element(gene1, "genie")[0].gene.label
        gene_3prime = self._get_gene_element(gene2, "genie")[0].gene.label

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(
                genie_row.get_column("Site1_Chromosome").item(), rb
            ),
            seg_end_genomic=int(genie_row.get_column("Site1_Position").item()),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(
                genie_row.get_column("Site2_Chromosome").item(), rb
            ),
            seg_start_genomic=int(genie_row.get_column("Site2_Position").item()),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(
            genie_row.get_column("Site1_Chromosome").item(),
            genie_row.get_column("Site2_Chromosome").item(),
            genie_row.get_column("Annotation").item(),
        )
        rf = bool(genie_row.get_column("Site2_Effect_on_Frame").item() == "in frame")
        return self._format_fusion(
            gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce, rf
        )
