"""Model for fusion class"""

from abc import ABC
from enum import Enum
from typing import Annotated, Any, Literal

from civicpy.civic import MolecularProfile
from cool_seq_tool.schemas import Strand
from ga4gh.core.models import MappableConcept
from ga4gh.vrs.models import (
    LiteralSequenceExpression,
    SequenceLocation,
)
from gene.schemas import CURIE_REGEX
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    StrictBool,
    StrictInt,
    StrictStr,
    StringConstraints,
    model_validator,
)

LINKER_REGEX = r"\|([atcg]+)\|"


class BaseModelForbidExtra(BaseModel, extra="forbid"):
    """Base Pydantic model class with extra values forbidden."""


class FUSORTypes(str, Enum):
    """Define FUSOR object type values."""

    FUNCTIONAL_DOMAIN = "FunctionalDomain"
    TRANSCRIPT_SEGMENT_ELEMENT = "TranscriptSegmentElement"
    TEMPLATED_SEQUENCE_ELEMENT = "TemplatedSequenceElement"
    LINKER_SEQUENCE_ELEMENT = "LinkerSequenceElement"
    GENE_ELEMENT = "GeneElement"
    UNKNOWN_GENE_ELEMENT = "UnknownGeneElement"
    MULTIPLE_POSSIBLE_GENES_ELEMENT = "MultiplePossibleGenesElement"
    BREAKPOINT_COVERAGE = "BreakpointCoverage"
    CONTIG_SEQUENCE = "ContigSequence"
    ANCHORED_READS = "AnchoredReads"
    SPLIT_READS = "SplitReads"
    SPANNING_READS = "SpanningReads"
    READ_DATA = "ReadData"
    REGULATORY_ELEMENT = "RegulatoryElement"
    CATEGORICAL_FUSION = "CategoricalFusion"
    ASSAYED_FUSION = "AssayedFusion"
    CAUSATIVE_EVENT = "CausativeEvent"


class AdditionalFields(str, Enum):
    """Define possible fields that can be added to Fusion object."""

    SEQUENCE_ID = "sequence_id"
    LOCATION_ID = "location_id"


class DomainStatus(str, Enum):
    """Define possible statuses of functional domains."""

    LOST = "lost"
    PRESERVED = "preserved"


class FunctionalDomain(BaseModel):
    """Define FunctionalDomain class"""

    type: Literal[FUSORTypes.FUNCTIONAL_DOMAIN] = FUSORTypes.FUNCTIONAL_DOMAIN
    status: DomainStatus
    associatedGene: MappableConcept
    id: Annotated[str, StringConstraints(pattern=CURIE_REGEX)] | None
    label: StrictStr | None = None
    sequenceLocation: SequenceLocation | None = None

    model_config = ConfigDict(
        populate_by_name=True,
        json_schema_extra={
            "example": {
                "type": "FunctionalDomain",
                "status": "lost",
                "label": "Tyrosine-protein kinase, catalytic domain",
                "id": "interpro:IPR020635",
                "associatedGene": {
                    "primaryCoding": {
                        "id": "hgnc:8031",
                        "code": "HGNC:8031",
                        "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                    },
                    "name": "NTRK1",
                    "conceptType": "Gene",
                },
                "sequenceLocation": {
                    "id": "ga4gh:SL.ywhUSfEUrwG0E29Q3c47bbuc6gkqTGlO",
                    "start": 510,
                    "end": 781,
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NC_000022.11",
                        "type": "SequenceReference",
                        "refgetAccession": "SQ.7B7SHsmchAR0dFcDCuSFjJAo7tX87krQ",
                    },
                },
            }
        },
    )


class StructuralElementType(str, Enum):
    """Define possible structural element type values."""

    TRANSCRIPT_SEGMENT_ELEMENT = FUSORTypes.TRANSCRIPT_SEGMENT_ELEMENT.value
    TEMPLATED_SEQUENCE_ELEMENT = FUSORTypes.TEMPLATED_SEQUENCE_ELEMENT.value
    LINKER_SEQUENCE_ELEMENT = FUSORTypes.LINKER_SEQUENCE_ELEMENT.value
    GENE_ELEMENT = FUSORTypes.GENE_ELEMENT.value
    UNKNOWN_GENE_ELEMENT = FUSORTypes.UNKNOWN_GENE_ELEMENT.value
    MULTIPLE_POSSIBLE_GENES_ELEMENT = FUSORTypes.MULTIPLE_POSSIBLE_GENES_ELEMENT.value


class BaseStructuralElement(ABC, BaseModel):
    """Define base structural element class."""

    type: StructuralElementType


class BreakpointCoverage(BaseStructuralElement):
    """Define BreakpointCoverage class.

    This class models breakpoint coverage, or the number of fragments
    that are retained near the breakpoint for a fusion partner
    """

    type: Literal[FUSORTypes.BREAKPOINT_COVERAGE] = FUSORTypes.BREAKPOINT_COVERAGE
    fragmentCoverage: int = Field(ge=0)

    model_config = ConfigDict(
        json_schema_extra={
            "example": {"type": "BreakpointCoverage", "fragmentCoverage": 180}
        }
    )


class ContigSequence(BaseStructuralElement):
    """Define ContigSequence class.

    This class models the assembled contig sequence that supports the reported fusion
    event
    """

    type: Literal[FUSORTypes.CONTIG_SEQUENCE] = FUSORTypes.CONTIG_SEQUENCE
    contig: Annotated[
        str,
        StringConstraints(
            strip_whitespace=True,
            to_upper=True,
            pattern=r"^(?:[^A-Za-z0-9]|[ACTGactg])*$",
        ),
    ]

    model_config = ConfigDict(
        json_schema_extra={
            "example": {"type": "ContigSequence", "contig": "GTACTACTGATCTAGCATCTAGTA"}
        }
    )


class AnchoredReads(BaseStructuralElement):
    """Define AnchoredReads class

    This class can be used to report the number of reads that span the
    fusion junction. This is used at the TranscriptSegment level, as it
    indicates the transcript where the longer segment of the read is found
    """

    type: Literal[FUSORTypes.ANCHORED_READS] = FUSORTypes.ANCHORED_READS
    reads: int = Field(ge=0)


class SplitReads(BaseStructuralElement):
    """Define SplitReads class.

    This class models the number of reads that cover the junction bewteen the
    detected partners in the fusion
    """

    type: Literal[FUSORTypes.SPLIT_READS] = FUSORTypes.SPLIT_READS
    splitReads: int = Field(ge=0)

    model_config = ConfigDict(
        json_schema_extra={"example": {"type": "SplitReads", "splitReads": 100}}
    )


class SpanningReads(BaseStructuralElement):
    """Define Spanning Reads class.

    This class models the number of pairs of reads that support the reported fusion
    event
    """

    type: Literal[FUSORTypes.SPANNING_READS] = FUSORTypes.SPANNING_READS
    spanningReads: int = Field(ge=0)

    model_config = ConfigDict(
        json_schema_extra={"example": {"type": "SpanningReads", "spanningReads": 100}}
    )


class ReadData(BaseStructuralElement):
    """Define ReadData class.

    This class is used at the AssayedFusion level when a fusion caller reports
    metadata describing sequencing reads for the fusion event
    """

    type: Literal[FUSORTypes.READ_DATA] = FUSORTypes.READ_DATA
    split: SplitReads | None = None
    spanning: SpanningReads | None = None

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "ReadData",
                "split": {"type": "SplitReads", "splitReads": 100},
                "spanning": {"type": "SpanningReads", "spanningReads": 80},
            }
        }
    )


class TranscriptSegmentElement(BaseStructuralElement):
    """Define TranscriptSegment class"""

    type: Literal[FUSORTypes.TRANSCRIPT_SEGMENT_ELEMENT] = (
        FUSORTypes.TRANSCRIPT_SEGMENT_ELEMENT
    )
    transcript: Annotated[str, StringConstraints(pattern=CURIE_REGEX)]
    strand: Strand
    exonStart: StrictInt | None = None
    exonStartOffset: StrictInt | None = 0
    exonEnd: StrictInt | None = None
    exonEndOffset: StrictInt | None = 0
    gene: MappableConcept
    elementGenomicStart: SequenceLocation | None = None
    elementGenomicEnd: SequenceLocation | None = None
    coverage: BreakpointCoverage | None = None
    anchoredReads: AnchoredReads | None = None

    @model_validator(mode="before")
    def check_exons(cls, values):
        """Check that at least one of {``exonStart``, ``exonEnd``} is set.
        If set, check that the corresponding ``elementGenomic`` field is set.
        If not set, set corresponding offset to ``None``

        """
        msg = "Must give values for either `exonStart`, `exonEnd`, or both"
        exon_start = values.get("exonStart")
        exon_end = values.get("exonEnd")
        if (exon_start is None) and (exon_end is None):
            raise ValueError(msg)

        if exon_start:
            if not values.get("elementGenomicStart"):
                msg = "Must give `elementGenomicStart` if `exonStart` is given"
                raise ValueError(msg)
        else:
            values["exonStartOffset"] = None

        if exon_end:
            if not values.get("elementGenomicEnd"):
                msg = "Must give `elementGenomicEnd` if `exonEnd` is given"
                raise ValueError(msg)
        else:
            values["exonEndOffset"] = None
        return values

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_152263.3",
                "strand": -1,
                "exonStart": 1,
                "exonStartOffset": 0,
                "exonEnd": 8,
                "exonEndOffset": 0,
                "gene": {
                    "primaryCoding": {
                        "id": "hgnc:12012",
                        "code": "HGNC:12012",
                        "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                    },
                    "conceptType": "Gene",
                    "name": "TPM3",
                },
                "elementGenomicStart": {
                    "id": "ga4gh:SL.Q8vkGp7_xR9vI0PQ7g1IvUUeQ4JlJG8l",
                    "digest": "Q8vkGp7_xR9vI0PQ7g1IvUUeQ4JlJG8l",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NC_000001.11",
                        "type": "SequenceReference",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                    },
                    "end": 154192135,
                },
                "elementGenomicEnd": {
                    "id": "ga4gh:SL.Lnne0bSsgjzmNkKsNnXg98FeJSrDJuLb",
                    "digest": "Lnne0bSsgjzmNkKsNnXg98FeJSrDJuLb",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NC_000001.11",
                        "type": "SequenceReference",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                    },
                    "start": 154170399,
                },
                "coverage": {
                    "type": "BreakpointCoverage",
                    "fragmentCoverage": 185,
                },
                "anchoredReads": {
                    "type": "AnchoredReads",
                    "reads": 100,
                },
            }
        },
    )


class LinkerElement(BaseStructuralElement, extra="forbid"):
    """Define Linker class (linker sequence)"""

    type: Literal[FUSORTypes.LINKER_SEQUENCE_ELEMENT] = (
        FUSORTypes.LINKER_SEQUENCE_ELEMENT
    )
    linkerSequence: LiteralSequenceExpression

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "LinkerSequenceElement",
                "linkerSequence": {
                    "id": "sequence:ACGT",
                    "type": "LiteralSequenceExpression",
                    "sequence": "ACGT",
                },
            }
        },
    )


class TemplatedSequenceElement(BaseStructuralElement):
    """Define Templated Sequence Element class.

    A templated sequence is a contiguous genomic sequence found in the gene
    product.
    """

    type: Literal[FUSORTypes.TEMPLATED_SEQUENCE_ELEMENT] = (
        FUSORTypes.TEMPLATED_SEQUENCE_ELEMENT
    )
    region: SequenceLocation
    strand: Strand

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "TemplatedSequenceElement",
                "region": {
                    "id": "ga4gh:SL.q_LeFVIakQtxnGHgxC4yehpLUxd6QsEr",
                    "type": "SequenceLocation",
                    "start": 44908821,
                    "end": 44908822,
                    "sequenceReference": {
                        "id": "refseq:NC_000012.12",
                        "refgetAccession": "SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl",
                    },
                },
                "strand": 1,
            }
        },
    )


class GeneElement(BaseStructuralElement):
    """Define Gene Element class."""

    type: Literal[FUSORTypes.GENE_ELEMENT] = FUSORTypes.GENE_ELEMENT
    gene: MappableConcept

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "GeneElement",
                "gene": {
                    "primaryCoding": {
                        "id": "hgnc:1097",
                        "code": "HGNC:1097",
                        "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                    },
                    "name": "BRAF",
                    "conceptType": "Gene",
                },
            }
        },
    )


class UnknownGeneElement(BaseStructuralElement):
    """Define UnknownGene class.

    This is primarily intended to represent a
    partner in the result of a fusion partner-agnostic assay, which identifies
    the absence of an expected gene. For example, a FISH break-apart probe may
    indicate rearrangement of an MLL gene, but by design, the test cannot
    provide the identity of the new partner. In this case, we would associate
    any clinical observations from this patient with the fusion of MLL with
    an UnknownGene element.
    """

    type: Literal[FUSORTypes.UNKNOWN_GENE_ELEMENT] = FUSORTypes.UNKNOWN_GENE_ELEMENT

    model_config = ConfigDict(
        json_schema_extra={"example": {"type": "UnknownGeneElement"}},
    )


class MultiplePossibleGenesElement(BaseStructuralElement):
    """Define MultiplePossibleGenesElement class.

    This is primarily intended to
    represent a partner in a categorical fusion, typifying generalizable
    characteristics of a class of fusions such as retained or lost regulatory elements
    and/or functional domains, often curated from biomedical literature for use in
    genomic knowledgebases. For example, EWSR1 rearrangements are often found in Ewing
    and Ewing-like small round cell sarcomas, regardless of the partner gene.
    We would associate this assertion with the fusion of EWSR1 with a
    MultiplePossibleGenesElement.
    """

    type: Literal[FUSORTypes.MULTIPLE_POSSIBLE_GENES_ELEMENT] = (
        FUSORTypes.MULTIPLE_POSSIBLE_GENES_ELEMENT
    )

    model_config = ConfigDict(
        json_schema_extra={"example": {"type": "MultiplePossibleGenesElement"}},
    )


class RegulatoryClass(str, Enum):
    """Define possible classes of Regulatory Elements.

    Options are the possible values for ``/regulatory_class`` value property in the
    `INSDC controlled vocabulary <https://www.insdc.org/controlled-vocabulary-regulatoryclass>`_.
    """

    ATTENUATOR = "attenuator"
    CAAT_SIGNAL = "caat_signal"
    ENHANCER = "enhancer"
    ENHANCER_BLOCKING_ELEMENT = "enhancer_blocking_element"
    GC_SIGNAL = "gc_signal"
    IMPRINTING_CONTROL_REGION = "imprinting_control_region"
    INSULATOR = "insulator"
    LOCUS_CONTROL_REGION = "locus_control_region"
    MINUS_35_SIGNAL = "minus_35_signal"
    MINUS_10_SIGNAL = "minus_10_signal"
    POLYA_SIGNAL_SEQUENCE = "polya_signal_sequence"
    PROMOTER = "promoter"
    RESPONSE_ELEMENT = "response_element"
    RIBOSOME_BINDING_SITE = "ribosome_binding_site"
    RIBOSWITCH = "riboswitch"
    SILENCER = "silencer"
    TATA_BOX = "tata_box"
    TERMINATOR = "terminator"
    OTHER = "other"


class RegulatoryElement(BaseModel):
    """Define RegulatoryElement class.

    ``featureId`` would ideally be constrained as a CURIE, but Encode, our preferred
    feature ID source, doesn't currently have a registered CURIE structure for ``EH_``
    identifiers. Consequently, we permit any kind of free text.
    """

    type: Literal[FUSORTypes.REGULATORY_ELEMENT] = FUSORTypes.REGULATORY_ELEMENT
    regulatoryClass: RegulatoryClass
    featureId: str | None = None
    associatedGene: MappableConcept | None = None
    featureLocation: SequenceLocation | None = None

    @model_validator(mode="before")
    def ensure_min_values(cls, values):
        """Ensure that one of {`featureId`, `featureLocation`}, and/or
        `associatedGene` is set.
        """
        if not (
            bool(values.get("featureId")) ^ bool(values.get("featureLocation"))
        ) and not (values.get("associatedGene")):
            msg = (
                "Must set 1 of {`featureId`, `associatedGene`} and/or `featureLocation`"
            )
            raise ValueError(msg)
        return values

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "RegulatoryElement",
                "regulatoryClass": "promoter",
                "featureLocation": {
                    "id": "ga4gh:SL.9hqdPDfXC-m_t_bDH75FZHfaM6OKDtRw",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NC_000001.11",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                    },
                    "start": 155593,
                    "end": 155610,
                },
            }
        },
    )


class FusionType(str, Enum):
    """Specify possible Fusion types."""

    CATEGORICAL_FUSION = FUSORTypes.CATEGORICAL_FUSION.value
    ASSAYED_FUSION = FUSORTypes.ASSAYED_FUSION.value

    @classmethod
    def values(cls) -> set:
        """Provide all possible enum values."""
        return {c.value for c in cls}


class AbstractFusion(BaseModel, ABC):
    """Define Fusion class"""

    type: FusionType
    regulatoryElement: RegulatoryElement | None = None
    structure: list[BaseStructuralElement]
    readingFramePreserved: StrictBool | None = None
    viccNomenclature: StrictStr | None = None

    @classmethod
    def _access_object_attr(
        cls,
        obj: dict | BaseModel,
        attr_name: str,
    ) -> Any | None:  # noqa: ANN401
        """Help enable safe access of object properties while performing validation for
        Pydantic class objects.

        Because the validator could be handling either
        existing Pydantic class objects, or candidate dictionaries, we need a flexible
        accessor function.

        :param obj: object to access
        :param attr_name: name of attribute to retrieve
        :return: attribute if successful retrieval, otherwise None
        :raise ValueError: if object doesn't have properties (ie it's not a dict or Pydantic
            model)
        """
        if isinstance(obj, BaseModel):
            try:
                return obj.__getattribute__(attr_name)
            except AttributeError:
                return None
        elif isinstance(obj, dict):
            return obj.get(attr_name)
        else:
            msg = "Unrecognized type, should only pass entities with properties"
            raise ValueError(msg)

    @classmethod
    def _fetch_gene_id(
        cls,
        obj: dict | BaseModel,
        alt_field: str | None = None,
    ) -> str | None:
        """Get gene ID if element includes a gene annotation.

        :param obj: element to fetch gene from. Might not contain a gene (e.g. it's a
            TemplatedSequenceElement) so we have to use safe checks to fetch.
        :param alt_field: the field to fetch the gene from, if it is not called "gene" (ex: associatedGene instead)
        :return: gene ID if gene is defined
        """
        gene_info = cls._access_object_attr(obj, alt_field if alt_field else "gene")
        if gene_info:
            gene_id = cls._access_object_attr(gene_info, "primaryCoding")
            if isinstance(gene_id, str):
                return gene_id
            gene_id = cls._access_object_attr(gene_id, "id")
            if gene_id:
                return gene_id
        return None

    @model_validator(mode="before")
    def enforce_abc(cls, values):
        """Ensure only subclasses can be instantiated."""
        if cls.__name__ == "AbstractFusion":
            msg = "Cannot instantiate Fusion abstract class"
            raise ValueError(msg)
        return values

    @model_validator(mode="before")
    def enforce_element_quantities(cls, values):
        """Ensure minimum # of elements, and require > 1 unique genes.

        To validate the unique genes rule, we extract gene IDs from the elements that
        designate genes, and take the number of total elements. If there is only one
        unique gene ID, and there are no non-gene-defining elements (such as
        an unknown partner), then we raise an error.
        """
        qt_error_msg = (
            "Fusions must contain >= 2 structural elements, or >=1 structural element "
            "and a regulatory element"
        )
        structure = values.get("structure", [])
        if not structure:
            raise ValueError(qt_error_msg)
        num_structure = len(structure)
        reg_element = values.get("regulatoryElement")
        if (num_structure + bool(reg_element)) < 2:
            raise ValueError(qt_error_msg)

        uq_gene_msg = "Fusions must form a chimeric transcript from two or more genes, or a novel interaction between a rearranged regulatory element with the expressed product of a partner gene."
        gene_ids = []
        if reg_element:
            gene_id = cls._fetch_gene_id(obj=reg_element, alt_field="associatedGene")
            if gene_id:
                gene_ids.append(gene_id)

        for element in structure:
            gene_id = cls._fetch_gene_id(obj=element)
            if gene_id:
                gene_ids.append(gene_id)

        unique_gene_ids = set(gene_ids)
        if len(unique_gene_ids) == 1 and len(gene_ids) == (
            num_structure + bool(reg_element)
        ):
            raise ValueError(uq_gene_msg)
        return values

    @model_validator(mode="after")
    def structure_ends(cls, values):
        """Ensure start/end elements are of legal types and have fields required by
        their position.
        """
        elements = values.structure
        error_messages = []
        if isinstance(elements[0], TranscriptSegmentElement):
            if elements[0].exonEnd is None and not values.regulatoryElement:
                msg = "5' TranscriptSegmentElement fusion partner must contain ending exon position"
                error_messages.append(msg)
        elif isinstance(elements[0], LinkerElement):
            msg = "First structural element cannot be LinkerSequence"
            error_messages.append(msg)

        if len(elements) > 2:
            for element in elements[1:-1]:
                if isinstance(element, TranscriptSegmentElement) and (
                    element.exonStart is None or element.exonEnd is None
                ):
                    msg = "Connective TranscriptSegmentElement must include both start and end positions"
                    error_messages.append(msg)
        if isinstance(elements[-1], TranscriptSegmentElement) and (
            elements[-1].exonStart is None
        ):
            msg = "3' fusion partner junction must include " "starting position"
            error_messages.append(msg)
        if error_messages:
            raise ValueError("\n".join(error_messages))
        return values


class Evidence(str, Enum):
    """Form of evidence supporting identification of the fusion."""

    OBSERVED = "observed"
    INFERRED = "inferred"


class Assay(BaseModelForbidExtra):
    """Information pertaining to the assay used in identifying the fusion."""

    type: Literal["Assay"] = "Assay"
    assayName: StrictStr | None = None
    assayId: Annotated[str, StringConstraints(pattern=CURIE_REGEX)] | None = None
    methodUri: Annotated[str, StringConstraints(pattern=CURIE_REGEX)] | None = None
    fusionDetection: Evidence | None = None

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "methodUri": "pmid:33576979",
                "assayId": "obi:OBI_0003094",
                "assayName": "fluorescence in-situ hybridization assay",
                "fusionDetection": "inferred",
            }
        }
    )


AssayedFusionElement = Annotated[
    TranscriptSegmentElement
    | GeneElement
    | TemplatedSequenceElement
    | LinkerElement
    | UnknownGeneElement
    | ContigSequence
    | ReadData,
    Field(discriminator="type"),
]


class EventType(str, Enum):
    """Permissible values for describing the underlying causative event driving an
    assayed fusion.
    """

    REARRANGEMENT = "rearrangement"
    READ_THROUGH = "read-through"
    TRANS_SPLICING = "trans-splicing"


class CausativeEvent(BaseModelForbidExtra):
    """Define causative event information for a fusion.

    The evaluation of a fusion may be influenced by the underlying mechanism that
    generated the fusion. Often this will be a DNA rearrangement, but it could also be
    a read-through or trans-splicing event.
    """

    type: Literal[FUSORTypes.CAUSATIVE_EVENT] = FUSORTypes.CAUSATIVE_EVENT
    eventType: EventType
    eventDescription: StrictStr | None = None

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "CausativeEvent",
                "eventType": "rearrangement",
                "eventDescription": "chr2:g.pter_8,247,756::chr11:g.15,825,273_cen_qter (der11) and chr11:g.pter_15,825,272::chr2:g.8,247,757_cen_qter (der2)",
            }
        },
    )


class AssayedFusion(AbstractFusion):
    """Assayed gene fusions from biological specimens are directly detected using
    RNA-based gene fusion assays, or alternatively may be inferred from genomic
    rearrangements detected by whole genome sequencing or by coarser-scale cytogenomic
    assays. Example: an EWSR1 fusion inferred from a breakapart FISH assay.
    """

    type: Literal[FUSORTypes.ASSAYED_FUSION] = FUSORTypes.ASSAYED_FUSION
    structure: list[AssayedFusionElement]
    causativeEvent: CausativeEvent | None = None
    assay: Assay | None = None
    contig: ContigSequence | None = None
    readData: ReadData | None = None

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "AssayedFusion",
                "causativeEvent": {
                    "type": "CausativeEvent",
                    "eventType": "rearrangement",
                    "eventDescription": "chr2:g.pter_8,247,756::chr11:g.15,825,273_cen_qter (der11) and chr11:g.pter_15,825,272::chr2:g.8,247,757_cen_qter (der2)",
                },
                "assay": {
                    "type": "Assay",
                    "methodUri": "pmid:33576979",
                    "assayId": "obi:OBI_0003094",
                    "assayName": "fluorescence in-situ hybridization assay",
                    "fusionDetection": "inferred",
                },
                "contig": {
                    "type": "ContigSequence",
                    "contig": "GTACTACTGATCTAGCATCTAGTA",
                },
                "readData": {
                    "type": "ReadData",
                    "split": {
                        "type": "SplitReads",
                        "splitReads": 100,
                    },
                    "spanning": {
                        "type": "SpanningReads",
                        "spanningReads": 80,
                    },
                },
                "structure": [
                    {
                        "type": "GeneElement",
                        "gene": {
                            "conceptType": "Gene",
                            "primaryCoding": {
                                "id": "hgnc:3058",
                                "code": "HGNC:3058",
                                "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                            },
                            "name": "EWSR1",
                        },
                    },
                    {"type": "UnknownGeneElement"},
                ],
            }
        },
    )


CategoricalFusionElement = Annotated[
    TranscriptSegmentElement
    | GeneElement
    | TemplatedSequenceElement
    | LinkerElement
    | MultiplePossibleGenesElement,
    Field(discriminator="type"),
]


class CategoricalFusion(AbstractFusion):
    """Categorical gene fusions are generalized concepts representing a class
    of fusions by their shared attributes, such as retained or lost regulatory
    elements and/or functional domains, and are typically curated from the
    biomedical literature for use in genomic knowledgebases.
    """

    type: Literal[FUSORTypes.CATEGORICAL_FUSION] = FUSORTypes.CATEGORICAL_FUSION
    criticalFunctionalDomains: list[FunctionalDomain] | None = None
    structure: list[CategoricalFusionElement]
    civicMolecularProfiles: list[MolecularProfile] | None = None

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        json_schema_extra={
            "example": {
                "type": "CategoricalFusion",
                "readingFramePreserved": True,
                "criticalFunctionalDomains": [
                    {
                        "type": "FunctionalDomain",
                        "status": "lost",
                        "label": "cystatin domain",
                        "id": "interpro:IPR000010",
                        "associatedGene": {
                            "primaryCoding": {
                                "id": "hgnc:2743",
                                "code": "HGNC:2743",
                                "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                            },
                            "name": "CST1",
                            "conceptType": "Gene",
                        },
                    }
                ],
                "structure": [
                    {
                        "type": "TranscriptSegmentElement",
                        "transcript": "refseq:NM_152263.3",
                        "strand": -1,
                        "exonStart": 1,
                        "exonStartOffset": 0,
                        "exonEnd": 8,
                        "exonEndOffset": 0,
                        "gene": {
                            "primaryCoding": {
                                "id": "hgnc:12012",
                                "code": "HGNC:12012",
                                "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                            },
                            "conceptType": "Gene",
                            "name": "TPM3",
                        },
                        "elementGenomicStart": {
                            "id": "ga4gh:SL.Q8vkGp7_xR9vI0PQ7g1IvUUeQ4JlJG8l",
                            "digest": "Q8vkGp7_xR9vI0PQ7g1IvUUeQ4JlJG8l",
                            "type": "SequenceLocation",
                            "sequenceReference": {
                                "id": "refseq:NC_000001.11",
                                "type": "SequenceReference",
                                "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                            },
                            "end": 154192135,
                        },
                        "elementGenomicEnd": {
                            "id": "ga4gh:SL.Lnne0bSsgjzmNkKsNnXg98FeJSrDJuLb",
                            "digest": "Lnne0bSsgjzmNkKsNnXg98FeJSrDJuLb",
                            "type": "SequenceLocation",
                            "sequenceReference": {
                                "id": "refseq:NC_000001.11",
                                "type": "SequenceReference",
                                "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                            },
                            "start": 154170399,
                        },
                    },
                    {
                        "type": "GeneElement",
                        "gene": {
                            "primaryCoding": {
                                "id": "hgnc:427",
                                "code": "HGNC:427",
                                "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                            },
                            "name": "ALK",
                            "conceptType": "Gene",
                        },
                    },
                ],
                "regulatoryElement": {
                    "type": "RegulatoryElement",
                    "regulatoryClass": "promoter",
                    "associatedGene": {
                        "conceptType": "Gene",
                        "primaryCoding": {
                            "id": "hgnc:1097",
                            "code": "HGNC:1097",
                            "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                        },
                        "name": "BRAF",
                    },
                },
            }
        },
    )


Fusion = CategoricalFusion | AssayedFusion
