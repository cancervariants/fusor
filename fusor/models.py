"""Model for fusion class"""
from typing import Optional, List, Union, Literal
from enum import Enum
from abc import ABC

from pydantic import BaseModel, validator, StrictInt, StrictBool, StrictStr, \
    Extra, ValidationError, root_validator
from ga4gh.vrsatile.pydantic import return_value
from ga4gh.vrsatile.pydantic.vrsatile_models import GeneDescriptor, \
    LocationDescriptor, SequenceDescriptor, CURIE
from ga4gh.vrsatile.pydantic.vrs_models import Sequence


class BaseModelForbidExtra(BaseModel):
    """Base model with extra fields forbidden."""

    class Config:
        """Configure class."""

        extra = Extra.forbid


class FUSORTypes(str, Enum):
    """Define FUSOR object type values."""

    FUNCTIONAL_DOMAIN = "FunctionalDomain"
    TRANSCRIPT_SEGMENT_COMPONENT = "TranscriptSegmentComponent"
    TEMPLATED_SEQUENCE_COMPONENT = "TemplatedSequenceComponent"
    LINKER_SEQUENCE_COMPONENT = "LinkerSequenceComponent"
    GENE_COMPONENT = "GeneComponent"
    UNKNOWN_GENE_COMPONENT = "UnknownGeneComponent"
    ANY_GENE_COMPONENT = "AnyGeneComponent"
    REGULATORY_ELEMENT = "RegulatoryElement"
    CATEGORICAL_FUSION = "CategoricalFusion"
    ASSAYED_FUSION = "AssayedFusion"


class AdditionalFields(str, Enum):
    """Define possible fields that can be added to Fusion object."""

    SEQUENCE_ID = "sequence_id"
    LOCATION_ID = "location_id"
    GENE_DESCRIPTOR = "gene_descriptor"


class DomainStatus(str, Enum):
    """Define possible statuses of functional domains."""

    LOST = "lost"
    PRESERVED = "preserved"


class FunctionalDomain(BaseModel):
    """Define FunctionalDomain class"""

    type: Literal[FUSORTypes.FUNCTIONAL_DOMAIN] = FUSORTypes.FUNCTIONAL_DOMAIN
    id: CURIE
    name: StrictStr
    status: DomainStatus
    gene_descriptor: GeneDescriptor
    location_descriptor: LocationDescriptor

    _get_id_val = validator("id", allow_reuse=True)(return_value)

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "FunctionalDomain",
                "status": "lost",
                "name": "Tyrosine-protein kinase, catalytic domain",
                "id": "interpro:IPR020635",
                "gene_descriptor": {
                    "id": "gene:NTRK1",
                    "gene_id": "hgnc:8031",
                    "label": "8031",
                    "type": "GeneDescriptor",
                },
                "location_descriptor": {
                    "id": "fusor.location_descriptor:NP_002520.2",
                    "type": "LocationDescriptor",
                    "location": {
                        "sequence_id": "ga4gh:SQ.vJvm06Wl5J7DXHynR9ksW7IK3_3jlFK6",  # noqa: E501
                        "type": "SequenceLocation",
                        "interval": {
                            "start": {
                                "type": "Number",
                                "value": 510
                            },
                            "end": {
                                "type": "Number",
                                "value": 781
                            }
                        }
                    }
                }
            }


class Component(ABC, BaseModel):
    """Define base component class."""


class TranscriptSegmentComponent(Component):
    """Define TranscriptSegment class"""

    type: Literal[FUSORTypes.TRANSCRIPT_SEGMENT_COMPONENT] = FUSORTypes.TRANSCRIPT_SEGMENT_COMPONENT  # noqa: E501
    transcript: CURIE
    exon_start: Optional[StrictInt]
    exon_start_offset: Optional[StrictInt] = 0
    exon_end: Optional[StrictInt]
    exon_end_offset: Optional[StrictInt] = 0
    gene_descriptor: GeneDescriptor
    component_genomic_start: Optional[LocationDescriptor]
    component_genomic_end: Optional[LocationDescriptor]

    @root_validator(pre=True)
    def check_exons(cls, values):
        """Check that at least one of {`exon_start`, `exon_end`} is set.
        If set, check that the corresponding `component_genomic` field is set.
        If not set, set corresponding offset to `None`

        """
        msg = "Must give values for either `exon_start`, `exon_end`, or both"
        exon_start = values.get("exon_start")
        exon_end = values.get("exon_end")
        assert exon_start or exon_end, msg

        if exon_start:
            msg = "Must give `component_genomic_start` if `exon_start` is given"  # noqa: E501
            assert values.get("component_genomic_start"), msg
        else:
            values["exon_start_offset"] = None

        if exon_end:
            msg = "Must give `component_genomic_end` if `exon_end` is given"
            assert values.get("component_genomic_end"), msg
        else:
            values["exon_end_offset"] = None
        return values

    _get_transcript_val = validator("transcript", allow_reuse=True)(return_value)  # noqa: E501

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "TranscriptSegmentComponent",
                "transcript": "refseq:NM_152263.3",
                "exon_start": 1,
                "exon_start_offset": 0,
                "exon_end": 8,
                "exon_end_offset": 0,
                "gene_descriptor": {
                    "type": "GeneDescriptor",
                    "id": "gene:TPM3",
                    "gene_id": "hgnc:12012",
                    "type": "GeneDescriptor",
                    "label": "TPM3",
                },
                "component_genomic_start": {
                    "id": "TPM3:exon1",
                    "type": "LocationDescriptor",
                    "location_id": "ga4gh:VSL.vyyyExx4enSZdWZr3z67-T8uVKH50uLi",  # noqa: E501
                    "location": {
                        "sequence_id": "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT",  # noqa: E501
                        "type": "SequenceLocation",
                        "interval": {
                            "start": {
                                "type": "Number",
                                "value": 154192135
                            },
                            "end": {
                                "type": "Number",
                                "value": 154192136
                            },
                            "type": "SequenceInterval"
                        }
                    }
                },
                "component_genomic_end": {
                    "id": "TPM3:exon8",
                    "type": "LocationDescriptor",
                    "location_id": "ga4gh:VSL._1bRdL4I6EtpBvVK5RUaXb0NN3k0gpqa",  # noqa: E501
                    "location": {
                        "sequence_id": "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT",  # noqa: E501
                        "type": "SequenceLocation",
                        "interval": {
                            "start": {
                                "type": "Number",
                                "value": 154170398
                            },
                            "end": {
                                "type": "Number",
                                "value": 154170399
                            },
                            "type": "SequenceInterval"
                        }
                    }
                }
            }


class LinkerComponent(Component):
    """Define Linker class (linker sequence)"""

    type: Literal[FUSORTypes.LINKER_SEQUENCE_COMPONENT] = FUSORTypes.LINKER_SEQUENCE_COMPONENT  # noqa: E501
    linker_sequence: SequenceDescriptor

    @validator("linker_sequence", pre=True)
    def validate_sequence(cls, v):
        """Enforce nucleotide base code requirements on sequence literals."""
        if isinstance(v, dict):
            try:
                v["sequence"] = v["sequence"].upper()
                seq = v["sequence"]
            except KeyError:
                raise TypeError
        elif isinstance(v, SequenceDescriptor):
            v.sequence = v.sequence.upper()
            seq = v.sequence
        else:
            raise TypeError

        try:
            Sequence(__root__=seq)
        except ValidationError:
            raise AssertionError("sequence does not match regex '^[A-Za-z*\\-]*$'")  # noqa: E501

        return v

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "LinkerSequenceComponent",
                "linker_sequence": {
                    "id": "sequence:ACGT",
                    "type": "SequenceDescriptor",
                    "sequence": "ACGT",
                    "residue_type": "SO:0000348"
                }
            }


class Strand(str, Enum):
    """Define possible values for strand"""

    POSITIVE = "+"
    NEGATIVE = "-"


class TemplatedSequenceComponent(BaseModel):
    """Define Templated Sequence Component class.
    A templated sequence is a contiguous genomic sequence found in the
    gene product
    """

    type: Literal[FUSORTypes.TEMPLATED_SEQUENCE_COMPONENT] = FUSORTypes.TEMPLATED_SEQUENCE_COMPONENT  # noqa: E501
    region: LocationDescriptor
    strand: Strand

    # TODO? add strand to sequencelocation, add chr property

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "TemplatedSequenceComponent",
                "region": {
                    "id": "chr12:44908821-44908822(+)",
                    "type": "LocationDescriptor",
                    "location_id": "ga4gh:VSL.AG54ZRBhg6pwpPLafF4KgaAHpdFio6l5",  # noqa: E501
                    "location": {
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl",  # noqa: E501
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 44908821},
                            "end": {"type": "Number", "value": 44908822}
                        },
                    },
                    "label": "chr12:44908821-44908822(+)"
                },
                "strand": "+"
            }


class GeneComponent(Component):
    """Define Gene component class."""

    type: Literal[FUSORTypes.GENE_COMPONENT] = FUSORTypes.GENE_COMPONENT
    gene_descriptor: GeneDescriptor

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "GeneComponent",
                "gene_descriptor": {
                    "id": "gene:BRAF",
                    "gene_id": "hgnc:1097",
                    "label": "BRAF",
                    "type": "GeneDescriptor",
                }
            }


class UnknownGeneComponent(BaseModel):
    """Define UnknownGene class. This is primarily intended to represent a
    partner in the result of a fusion partner-agnostic assay, which identifies
    the absence of an expected gene. For example, a FISH break-apart probe may
    indicate rearrangement of an MLL gene, but by design, the test cannot
    provide the identity of the new partner. In this case, we would associate
    any clinical observations from this patient with the fusion of MLL with
    an UnknownGene component.
    """

    type: Literal[FUSORTypes.UNKNOWN_GENE_COMPONENT] = FUSORTypes.UNKNOWN_GENE_COMPONENT  # noqa: E501

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "UnknownGeneComponent"
            }


class AnyGeneComponent(BaseModel):
    """Define AnyGene class. This is primarily intended to represent a partner
    in a categorical fusion, typifying generalizable characteristics of a class
    of fusions such as retained or lost regulatory elements and/or functional
    domains, often curated from biomedical literature for use in genomic
    knowledgebases. For example, EWSR1 rearrangements are often found in Ewing
    and Ewing-like small round cell sarcomas, regardless of the partner gene.
    We would associate this assertion with the fusion of EWSR1 with an
    AnyGene component.
    """

    type: Literal[FUSORTypes.ANY_GENE_COMPONENT] = FUSORTypes.ANY_GENE_COMPONENT  # noqa: E501

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "AnyGeneComponent"
            }


class Event(str, Enum):
    """Define Event class (causative event)"""

    REARRANGEMENT = "rearrangement"
    READTHROUGH = "read-through"
    TRANSSPLICING = "trans-splicing"


class RegulatoryElementType(str, Enum):
    """Define possible types of Regulatory Elements. Options are the possible
    values for /regulatory_class value property in the INSDC controlled
    vocabulary.
    https://www.insdc.org/controlled-vocabulary-regulatoryclass
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
    """Define RegulatoryElement class"""

    type: Literal[FUSORTypes.REGULATORY_ELEMENT] = FUSORTypes.REGULATORY_ELEMENT  # noqa: E501
    element_type: RegulatoryElementType
    element_reference: Optional[CURIE] = None
    associated_gene: Optional[GeneDescriptor] = None
    genomic_location: Optional[LocationDescriptor] = None

    _get_ref_id_val = validator("element_reference", allow_reuse=True)(return_value)  # noqa: E501

    @root_validator(pre=True)
    def ensure_min_values(cls, values):
        """Ensure one of {`element_reference`, `associated_gene`,
        `genomic_location`} is set.
        """
        if not any([values.get("element_reference"),
                    values.get("associated_gene"),
                    values.get("genomic_location")]):  # noqa: E501
            raise ValueError
        return values

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "RegulatoryElement",
                "element_type": "Promoter",
                "genomic_location": {
                    "type": "LocationDescriptor",
                    "location": {
                        "sequence_id": "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",  # noqa: E501
                        "type": "SequenceLocation",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {
                                "type": "Number",
                                "value": 155593
                            },
                            "end": {
                                "type": "Number",
                                "value": 155610
                            }
                        }
                    }
                }
            }


class FusionType(str, Enum):
    """Specify possible Fusion types."""

    CATEGORICAL_FUSION = FUSORTypes.CATEGORICAL_FUSION
    ASSAYED_FUSION = FUSORTypes.ASSAYED_FUSION


class Fusion(BaseModel, ABC):
    """Define Fusion class"""

    type: FusionType
    regulatory_elements: Optional[List[RegulatoryElement]]
    structural_components: List[Component]

    @root_validator(pre=True)
    def enforce_abc(cls, values):
        """Ensure only subclasses can be instantiated."""
        if cls.__name__ == "Fusion":
            raise ValueError("Cannot instantiate Fusion abstract class")
        return values

    @root_validator(skip_on_failure=True)
    def check_components_length(cls, values):
        """Ensure >=2 structural components + regulatory elements"""
        components = len(values.get("structural_components", []))
        if values.get("regulatory_elements"):
            elements = len(values["regulatory_elements"])
        else:
            elements = 0
        if (components < 1) or (components + elements < 2):
            raise ValueError("Provided fusion contains an insufficient number "
                             "of structural components and regulatory "
                             "elements.")
        else:
            return values

    @root_validator(skip_on_failure=True)
    def structural_components_ends(cls, values):
        """Ensure start/end components are of legal types and have fields
        required by their position.
        """
        components = values.get("structural_components")
        if isinstance(components[0], TranscriptSegmentComponent):
            if components[0].exon_end is None and not values["regulatory_elements"]:  # noqa: E501
                raise ValueError(
                    "5' TranscriptSegmentComponent fusion partner must "
                    "contain ending exon position"
                )
        elif isinstance(components[0], LinkerComponent):
            raise ValueError(
                "First structural component cannot be LinkerSequence")

        if len(components) > 2:
            for component in components[1:-1]:
                if isinstance(component, TranscriptSegmentComponent):
                    if component.exon_start is None or component.exon_end is None:  # noqa: E501
                        raise ValueError(
                            "Connective TranscriptSegmentComponents must "
                            "include both start and end positions"
                        )
        if isinstance(components[-1], TranscriptSegmentComponent):
            if components[-1].exon_start is None:
                raise ValueError("3' fusion partner junction must include "
                                 "starting position")
        return values


class Evidence(str, Enum):
    """Form of evidence supporting identification of the fusion."""

    OBSERVED = "observed"
    INFERRED = "inferred"


class MolecularAssay(BaseModelForbidExtra):
    """Information pertaining to the assay used in identifying the fusion."""

    assay_name: Optional[StrictStr]
    eco_id: Optional[CURIE]
    eco_label: Optional[StrictStr]
    method_uri: Optional[CURIE]

    _get_eco_id_val = validator("eco_id", allow_reuse=True)(return_value)
    _get_method_id_val = validator("method_uri", allow_reuse=True)(return_value)  # noqa: E501

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            schema["example"] = {
                "method_uri": "pmid:34974290",
                "eco_id": "ECO:0005629",
                "eco_label": "DNA affinity chromatography evidence used in manual assertion",  # noqa: E501
                "assay_name": "Breakapart FISH probe"
            }


AssayedFusionComponents = List[
    Union[
        TranscriptSegmentComponent, GeneComponent, TemplatedSequenceComponent,
        LinkerComponent, UnknownGeneComponent
    ]
]


class EventDescription(BaseModelForbidExtra):
    """Describe causative event."""

    pass  # TODO ???


class AssayedFusion(Fusion):
    """Assayed gene fusions from biological specimens are directly detected
    using RNA-based gene fusion assays, or alternatively may be inferred from
    genomic rearrangements detected by whole genome sequencing or by
    coarser-scale cytogenomic assays. Example: an EWSR1 fusion inferred
    from a breakapart FISH assay.
    """

    type: Literal[FUSORTypes.ASSAYED_FUSION] = FUSORTypes.ASSAYED_FUSION
    causative_event: Optional[Event]
    event_description: Optional[EventDescription]
    fusion_evidence: Optional[Evidence]
    molecular_assay: Optional[MolecularAssay]
    structural_components: AssayedFusionComponents

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "AssayedFusion",
                "causative_event": "rearrangement",
                "fusion_evidence": "observed",
                "molecular_assay": '{"ECO_id": "ECO:0007159"}',
                "structural_components": [
                    {
                        "component_type": "gene",
                        "gene_descriptor": {
                            "id": "gene:EWSR1",
                            "gene_id": "hgnc:3058",
                            "label": "EWSR1",
                            "type": "GeneDescriptor",
                        }
                    },
                    {
                        "component_type": "UnknownGene"
                    }
                ]
            }


CategoricalFusionComponents = List[Union[TranscriptSegmentComponent,
                                         GeneComponent,
                                         TemplatedSequenceComponent,
                                         LinkerComponent,
                                         AnyGeneComponent]]


class CategoricalFusion(Fusion):
    """Categorical gene fusions are generalized concepts representing a class
    of fusions by their shared attributes, such as retained or lost regulatory
    elements and/or functional domains, and are typically curated from the
    biomedical literature for use in genomic knowledgebases.
    """

    type: Literal[FUSORTypes.CATEGORICAL_FUSION] = FUSORTypes.CATEGORICAL_FUSION  # noqa: E501
    r_frame_preserved: Optional[StrictBool]
    functional_domains: Optional[List[FunctionalDomain]]
    structural_components: CategoricalFusionComponents

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "CategoricalFusion",
                "r_frame_preserved": True,
                "functional_domains": [
                    {
                        "status": "lost",
                        "name": "cystatin domain",
                        "id": "interpro:IPR000010",
                        "gene": {
                            "id": "gene:CST1",
                            "gene_id": "hgnc:2743",
                            "label": "CST1",
                            "type": "GeneDescriptor",
                        }
                    }
                ],
                "structural_components": [
                    {
                        "component_type": "TranscriptSegment",
                        "transcript": "refseq:NM_152263.3",
                        "exon_start": 1,
                        "exon_start_offset": 0,
                        "exon_end": 8,
                        "exon_end_offset": 0,
                        "gene": {
                            "id": "gene:TPM3",
                            "gene_id": "hgnc:12012",
                            "type": "GeneDescriptor",
                            "label": "TPM3",
                        },
                        "component_genomic_start": {
                            "id": "TPM3:exon1",
                            "type": "LocationDescriptor",
                            "location_id": "ga4gh:VSL.vyyyExx4enSZdWZr3z67-T8uVKH50uLi",  # noqa: E501
                            "location": {
                                "sequence_id": "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT",  # noqa: E501
                                "type": "SequenceLocation",
                                "interval": {
                                    "start": {
                                        "type": "Number",
                                        "value": 154192135
                                    },
                                    "end": {
                                        "type": "Number",
                                        "value": 154192136
                                    },
                                    "type": "SequenceInterval"
                                }
                            }
                        },
                        "component_genomic_end": {
                            "id": "TPM3:exon8",
                            "type": "LocationDescriptor",
                            "location_id": "ga4gh:VSL._1bRdL4I6EtpBvVK5RUaXb0NN3k0gpqa",  # noqa: E501
                            "location": {
                                "sequence_id": "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT",  # noqa: E501
                                "type": "SequenceLocation",
                                "interval": {
                                    "start": {
                                        "type": "Number",
                                        "value": 154170398
                                    },
                                    "end": {
                                        "type": "Number",
                                        "value": 154170399
                                    },
                                    "type": "SequenceInterval"
                                }
                            }
                        }
                    },
                    {
                        "component_type": "gene",
                        "gene": {
                            "id": "gene:ALK",
                            "type": "GeneDescriptor",
                            "gene_id": "hgnc:427",
                            "label": "ALK"
                        }
                    }
                ],
                "regulatory_elements": [
                    {
                        "element_type": "promoter",
                        "gene": {
                            "id": "gene:BRAF",
                            "type": "GeneDescriptor",
                            "gene_id": "hgnc:1097",
                            "label": "BRAF"
                        }
                    }
                ]
            }
