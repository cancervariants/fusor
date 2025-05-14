"""Schemas for outputs provided by different fusion callers"""

from abc import ABC
from enum import Enum
from typing import Literal

from civicpy.civic import ExonCoordinate, MolecularProfile
from pydantic import BaseModel, ConfigDict, Field


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


class KnowledgebaseList(str, Enum):
    """Define supported knowledgebases"""

    CIVIC = "CIVIC"


class FusionCaller(ABC, BaseModel):
    """ABC for fusion callers"""

    type: Caller
    model_config = ConfigDict(extra="allow")


class FusionKnowledgebase(ABC, BaseModel):
    """ABC for Fusion Knowledgebases"""

    type: KnowledgebaseList
    model_config = ConfigDict(extra="allow")


class JAFFA(FusionCaller):
    """Define parameters for JAFFA model"""

    type: Literal[Caller.JAFFA] = Caller.JAFFA
    fusion_genes: str = Field(
        ..., description="A string containing the two fusion partners"
    )
    chrom1: str = Field(
        ..., description="The chromosome indicated in the chrom1 column"
    )
    base1: int = Field(
        ..., description="The genomic position indicated in the base1 column"
    )
    chrom2: str = Field(
        ..., description="The chromosome indicated in the chrom2 column"
    )
    base2: int = Field(
        ..., description="The genomic position indicated in the base2 column"
    )
    rearrangement: bool = Field(
        ..., description=" A boolean indicating if a rearrangement occurred"
    )
    classification: str = Field(
        ..., description="The classification associated with the called fusion"
    )
    inframe: bool | str = Field(
        ...,
        description="A boolean or string indicating if the fusion occurred in-frame",
    )
    spanning_reads: int = Field(
        ...,
        description="The number of detected reads that span the junction between the two transcript. Although described as spanning reads, this aligns with our definition of split reads i.e. reads that have sequence belonging to the two fusion partners",
    )
    spanning_pairs: int = Field(
        ...,
        description="The number of detected reads that align entirely on either side of the breakpoint",
    )


class STARFusion(FusionCaller):
    """Define parameters for STAR-Fusion model"""

    type: Literal[Caller.STAR_FUSION] = Caller.STAR_FUSION
    left_gene: str = Field(..., description="The gene indicated in the LeftGene column")
    right_gene: str = Field(
        ..., description="The gene indicated in the RightGene column"
    )
    left_breakpoint: str = Field(
        ..., description="The gene indicated in the LeftBreakpoint column"
    )
    right_breakpoint: str = Field(
        ..., description="The gene indicated in the RightBreakpoint column"
    )
    annots: str = Field(..., description="The annotations associated with the fusion")
    junction_read_count: int = Field(
        ...,
        description="The number of RNA-seq fragments that split the junction between the two transcript segments (from STAR-Fusion documentation)",
    )
    spanning_frag_count: int = Field(
        ...,
        description="The number of RNA-seq fragments that encompass the fusion junction such that one read of the pair aligns to a different gene than the other paired-end read of that fragment (from STAR-Fusion documentation)",
    )


class FusionCatcher(FusionCaller):
    """Define parameters for FusionCatcher model"""

    type: Literal[Caller.FUSION_CATCHER] = Caller.FUSION_CATCHER
    five_prime_partner: str = Field(
        ..., description="Gene symbol for the 5' fusion partner"
    )
    three_prime_partner: str = Field(
        ..., description="Gene symbol for the 3' fusion partner"
    )
    five_prime_fusion_point: str = Field(
        ...,
        description="Chromosomal position for the 5' end of the fusion junction. This coordinate is 1-based",
    )
    three_prime_fusion_point: str = Field(
        ...,
        description="Chromosomal position for the 3' end of the fusion junction. This coordinate is 1-based",
    )
    predicted_effect: str = Field(
        ...,
        description="The predicted effect of the fusion event, created using annotation from the Ensembl database",
    )
    spanning_unique_reads: int = Field(
        ..., description="The number of unique reads that map on the fusion junction"
    )
    spanning_reads: int = Field(
        ..., description="The number of paired reads that support the fusion"
    )
    fusion_sequence: str = Field(
        ..., description="The inferred sequence around the fusion junction"
    )


class Arriba(FusionCaller):
    """Define parameters for Arriba model"""

    type: Literal[Caller.ARRIBA] = Caller.ARRIBA
    gene1: str = Field(..., description="The 5' gene fusion partner")
    gene2: str = Field(..., description="The 3' gene fusion partner")
    strand1: str = Field(
        ..., description="The strand information for the 5' gene fusion partner"
    )
    strand2: str = Field(
        ..., description="The strand information for the 3' gene fusion partner"
    )
    breakpoint1: str = Field(..., description="The chromosome and breakpoint for gene1")
    breakpoint2: str = Field(..., description="The chromosome and breakpoint for gene2")
    event_type: str = Field(
        ..., description=" An inference about the type of fusion event"
    )
    confidence: str = Field(
        ..., description="A metric describing the confidence of the fusion prediction"
    )
    direction1: str = Field(
        ...,
        description="A description that indicates if the transcript segment starts or ends at breakpoint1",
    )
    direction2: str = Field(
        ...,
        description="A description that indicates if the transcript segment starts or ends at breakpoint2",
    )
    rf: str = Field(
        ...,
        description="A description if the reading frame is preserved for the fusion",
    )
    split_reads1: int = Field(
        ..., description="Number of supporting split fragments with anchor in gene1"
    )
    split_reads2: int = Field(
        ..., description="Number of supporting split fragments with anchor in gene2"
    )
    discordant_mates: int = Field(
        ..., description="Number of discordant mates supporting the fusion"
    )
    coverage1: int = Field(
        ..., description="Number of fragments retained near breakpoint1"
    )
    coverage2: int = Field(
        ..., description="Number of fragments retained near breakpoint2"
    )
    fusion_transcript: str = Field(..., description="The assembled fusion transcript")


class Cicero(FusionCaller):
    """Define parameters for CICERO model"""

    type: Literal[Caller.CICERO] = Caller.CICERO
    gene_5prime: str = Field(..., description="The gene symbol for the 5' partner")
    gene_3prime: str = Field(..., description="The gene symbol for the 3' partner")
    chr_5prime: str = Field(..., description="The chromosome for the 5' partner")
    chr_3prime: str = Field(..., description="The chromosome for the 3' partner")
    pos_5prime: int = Field(
        ..., description="The genomic breakpoint for the 5' partner"
    )
    pos_3prime: int = Field(
        ..., description="The genomic breakpoint for the 3' partner"
    )
    sv_ort: str = Field(
        ...,
        description="Whether the mapping orientation of assembled contig (driven by structural variation) has confident biological meaning",
    )
    event_type: str = Field(
        ...,
        description="The structural variation event that created the called fusion",
    )
    reads_5prime: int = Field(
        ...,
        description="The number of reads that support the breakpoint for the 5' partner",
    )
    reads_3prime: int = Field(
        ...,
        description="The number of reads that support the breakpoint for the 3' partner",
    )
    coverage_5prime: int = Field(
        ..., description="The fragment coverage at the 5' breakpoint"
    )
    coverage_3prime: int = Field(
        ..., description="The fragment coverage at the 3' breakpoint"
    )
    contig: str = Field(..., description="The assembled contig sequence for the fusion")


class EnFusion(FusionCaller):
    """Define parameters for EnFusion model"""

    type: Literal[Caller.ENFUSION] = Caller.ENFUSION
    gene_5prime: str = Field(..., description="The 5' gene fusion partner")
    gene_3prime: str = Field(..., description="The 3' gene fusion partner")
    chr_5prime: int | str = Field(
        ..., description="The 5' gene fusion partner chromosome"
    )
    chr_3prime: int | str = Field(
        ..., description="The 3' gene fusion partner chromosome"
    )
    break_5prime: int = Field(
        ..., description="The 5' gene fusion partner genomic breakpoint"
    )
    break_3prime: int = Field(
        ..., description="The 3' gene fusion partner genomic breakpoint"
    )
    fusion_junction_sequence: str | None = Field(
        None, description="The sequence near the fusion junction"
    )


class Genie(FusionCaller):
    """Define parameters for Genie model"""

    type: Literal[Caller.GENIE] = Caller.GENIE
    site1_hugo: str = Field(..., description="The HUGO symbol reported at site 1")
    site2_hugo: str = Field(..., description="The HUGO symbol reported at site 2")
    site1_chrom: int = Field(..., description="The chromosome reported at site 1")
    site2_chrom: int = Field(..., description="The chromosome reported at site 2")
    site1_pos: int = Field(..., description="The breakpoint reported at site 1")
    site2_pos: int = Field(..., description="The breakpoint reported at site 2")
    annot: str = Field(..., description="The annotation for the fusion event")
    reading_frame: str = Field(
        ..., description="The reading frame status of the fusion"
    )


class CIVIC(FusionKnowledgebase):
    """Define parameters for CIVIC model"""

    type: Literal[KnowledgebaseList.CIVIC] = KnowledgebaseList.CIVIC
    model_config = ConfigDict(arbitrary_types_allowed=True)
    vicc_compliant_name: str = Field(
        ..., description="The VICC compliant name for the fusion"
    )
    five_prime_end_exon_coords: ExonCoordinate | None = Field(
        ..., description="Data for the end exon of 5' fusion partner"
    )
    three_prime_start_exon_coords: ExonCoordinate | None = Field(
        ..., description="Data for the start exon 3' fusion partner"
    )
    molecular_profiles: list[MolecularProfile] | None = Field(
        ..., description="The molecular profiles associated with the fusion"
    )
