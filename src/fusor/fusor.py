"""Module for modifying fusion objects."""

import logging
import re

from bioutils.accessions import coerce_namespace
from cool_seq_tool.app import CoolSeqTool
from cool_seq_tool.schemas import ResidueMode
from ga4gh.core import ga4gh_identify
from ga4gh.core.domain_models import Gene
from ga4gh.vrs import models
from ga4gh.vrs.models import LiteralSequenceExpression, SequenceLocation, SequenceString
from gene.database import AbstractDatabase as GeneDatabase
from gene.database import create_db
from gene.query import QueryHandler
from gene.schemas import CURIE
from pydantic import ValidationError

from fusor.exceptions import FUSORParametersException, IDTranslationException
from fusor.models import (
    AdditionalFields,
    Assay,
    AssayedFusion,
    AssayedFusionElements,
    BaseStructuralElement,
    CategoricalFusion,
    CategoricalFusionElements,
    CausativeEvent,
    DomainStatus,
    Evidence,
    FunctionalDomain,
    Fusion,
    FusionType,
    GeneElement,
    LinkerElement,
    MultiplePossibleGenesElement,
    RegulatoryClass,
    RegulatoryElement,
    Strand,
    StructuralElementType,
    TemplatedSequenceElement,
    TranscriptSegmentElement,
    UnknownGeneElement,
)
from fusor.nomenclature import (
    gene_nomenclature,
    reg_element_nomenclature,
    templated_seq_nomenclature,
    tx_segment_nomenclature,
)
from fusor.tools import translate_identifier

_logger = logging.getLogger(__name__)


class FUSOR:
    """Class for modifying fusion objects."""

    def __init__(
        self,
        gene_database: GeneDatabase | None = None,
        cool_seq_tool: CoolSeqTool | None = None,
    ) -> None:
        """Initialize FUSOR class.

        :param cool_seq_tool: Cool-Seq-Tool instance
        :param gene_database: gene normalizer database instance
        """
        if not gene_database:
            gene_database = create_db()
        self.gene_normalizer = QueryHandler(gene_database)

        if not cool_seq_tool:
            cool_seq_tool = CoolSeqTool()
        self.cool_seq_tool = cool_seq_tool
        self.seqrepo = self.cool_seq_tool.seqrepo_access.sr

    @staticmethod
    def _contains_element_type(kwargs: dict, elm_type: StructuralElementType) -> bool:
        """Check if fusion contains element of a specific type. Helper method for
        inferring fusion type.

        :param kwargs: keyword args given to fusion method
        :param elm_type: element type to match
        :return: True if at least one element of given type is found, False otherwise.
        """
        for c in kwargs["structure"]:
            if (isinstance(c, dict) and c.get("type") == elm_type) or (
                isinstance(c, BaseStructuralElement) and c.type == elm_type
            ):
                return True
        return False

    def fusion(self, fusion_type: FusionType | None = None, **kwargs) -> Fusion:
        """Construct fusion object.

        :param fusion_type: explicitly specify fusion type. Unnecessary if providing
            fusion object in keyword args that includes ``type`` attribute.
        :return: constructed fusion object if successful
        :raise: FUSORParametersException if fusion type unable to be determined,
            or if incorrect fusion parameters are provided
        """
        # try explicit type param
        explicit_type = kwargs.get("type")
        if not fusion_type and explicit_type:
            if explicit_type in FusionType.values():
                fusion_type = explicit_type
                kwargs.pop("type")
            else:
                msg = f"Invalid type parameter: {explicit_type}"
                raise FUSORParametersException(msg)
        fusion_fn = None
        if fusion_type:
            if fusion_type == FusionType.CATEGORICAL_FUSION:
                fusion_fn = self.categorical_fusion
            elif fusion_type == FusionType.ASSAYED_FUSION:
                fusion_fn = self.assayed_fusion
            else:
                msg = f"Invalid fusion_type parameter: {fusion_type}"
                raise FUSORParametersException(msg)
        else:
            # try to infer from provided attributes
            categorical_attributes = any(
                [
                    "criticalFunctionalDomains" in kwargs,
                    self._contains_element_type(
                        kwargs, StructuralElementType.MULTIPLE_POSSIBLE_GENES_ELEMENT
                    ),
                ]
            )
            assayed_attributes = any(
                [
                    "causativeEvent" in kwargs,
                    "assay" in kwargs,
                    self._contains_element_type(
                        kwargs, StructuralElementType.UNKNOWN_GENE_ELEMENT
                    ),
                ]
            )
            if categorical_attributes and assayed_attributes:
                msg = "Received conflicting attributes"
                raise FUSORParametersException(msg)
            if categorical_attributes and not assayed_attributes:
                fusion_fn = self.categorical_fusion
            elif assayed_attributes and not categorical_attributes:
                fusion_fn = self.assayed_fusion
        if fusion_fn is None:
            msg = "Unable to determine fusion type"
            raise FUSORParametersException(msg)
        try:
            return fusion_fn(**kwargs)
        except TypeError as e:
            msg = f"Unable to construct fusion with provided args: {e}"
            raise FUSORParametersException(msg) from e

    @staticmethod
    def categorical_fusion(
        structure: CategoricalFusionElements,
        regulatory_element: RegulatoryElement | None = None,
        critical_functional_domains: list[FunctionalDomain] | None = None,
        reading_frame_preserved: bool | None = None,
    ) -> CategoricalFusion:
        """Construct a categorical fusion object
        :param structure: elements constituting the fusion
        :param regulatory_element: affected regulatory element
        :param critical_functional_domains: lost or preserved functional domains
        :param reading_frame_preserved: ``True`` if reading frame is preserved.
            ``False`` otherwise
        :return: CategoricalFusion if construction successful
        :raise: FUSORParametersException if given incorrect fusion properties
        """
        try:
            fusion = CategoricalFusion(
                structure=structure,
                criticalFunctionalDomains=critical_functional_domains,
                readingFramePreserved=reading_frame_preserved,
                regulatoryElement=regulatory_element,
            )
        except ValidationError as e:
            raise FUSORParametersException(str(e)) from e
        return fusion

    @staticmethod
    def assayed_fusion(
        structure: AssayedFusionElements,
        causative_event: CausativeEvent | None = None,
        assay: Assay | None = None,
        regulatory_element: RegulatoryElement | None = None,
        reading_frame_preserved: bool | None = None,
    ) -> AssayedFusion:
        """Construct an assayed fusion object
        :param structure: elements constituting the fusion
        :param causativeEvent: event causing the fusion
        :param assay: how knowledge of the fusion was obtained
        :param regulatory_element: affected regulatory elements
        :param reading_frame_preserved: ``True`` if reading frame is preserved.
            ``False`` otherwise
        :return: Tuple containing optional AssayedFusion if construction successful,
            and any relevant validation warnings
        """
        try:
            fusion = AssayedFusion(
                structure=structure,
                regulatoryElement=regulatory_element,
                causativeEvent=causative_event,
                assay=assay,
                reading_frame_preserved=reading_frame_preserved,
            )
        except ValidationError as e:
            raise FUSORParametersException(str(e)) from e
        return fusion

    async def transcript_segment_element(
        self,
        tx_to_genomic_coords: bool = True,
        use_minimal_gene_descr: bool = True,
        seq_id_target_namespace: str | None = None,
        **kwargs,
    ) -> tuple[TranscriptSegmentElement | None, list[str] | None]:
        """Create transcript segment element

        :param tx_to_genomic_coords: `True` if going from transcript to genomic
            coordinates. ``False`` if going from genomic to transcript exon coordinates.
        :param use_minimal_gene_descr: `True` if minimal gene object
            (``id``, ``label``) will be used. ``False`` if
            gene-normalizer's entire gene object will be used
        :param seq_id_target_namespace: If want to use digest for ``sequence_id``, set
            this to the namespace you want the digest for. Otherwise, leave as ``None``.
        :param kwargs:
            If `tx_to_genomic_coords`, possible key word arguments:
                (From cool_seq_tool.transcript_to_genomic_coords)
                gene: Optional[str] = None, transcript: str = None,
                exon_start: Optional[int] = None,
                exon_start_offset: Optional[int] = 0,
                exon_end: Optional[int] = None,
                exon_end_offset: Optional[int] = 0
            else:
                (From cool_seq_tool.genomic_to_transcript_exon_coordinates)
                chromosome: Union[str, int], start: Optional[int] = None,
                end: Optional[int] = None, strand: Optional[int] = None,
                transcript: Optional[str] = None, gene: Optional[str] = None,
                residue_mode: ResidueMode = ResidueMode.RESIDUE
        :return: Transcript Segment Element, warning
        """
        if tx_to_genomic_coords:
            data = await self.cool_seq_tool.ex_g_coords_mapper.transcript_to_genomic_coordinates(
                **kwargs
            )
        else:
            if "chromosome" in kwargs and kwargs.get("chromosome") is None:
                msg = (
                    "`chromosome` is required when going from genomic to"
                    " transcript exon coordinates"
                )
                _logger.warning(msg)
                return None, [msg]
            residue_mode = kwargs.get("residue_mode")
            # TODO: Remove once fixed in cool_seq_tool - need to verify that this code isn't still needed
            if residue_mode != ResidueMode.INTER_RESIDUE:
                start = kwargs.get("start")
                kwargs["start"] = start - 1 if start is not None else None
                kwargs["residue_mode"] = "inter-residue"
            chromosome = kwargs.get("chromosome")
            # if chromosome is a string, assume it's an accession, fix it for the kwargs since CST expects
            if type(chromosome) is str:
                kwargs["alt_ac"] = chromosome
            data = await self.cool_seq_tool.ex_g_coords_mapper.genomic_to_transcript_exon_coordinates(
                **kwargs
            )

        if data.genomic_data is None:
            return None, data.warnings

        genomic_data = data.genomic_data
        # should we be doing this before the call to cool_seq_tool? if the namespace prefix is present, that call will fail
        genomic_data.transcript = coerce_namespace(genomic_data.transcript)

        normalized_gene_response = self._normalized_gene(
            genomic_data.gene, use_minimal_gene_descr=use_minimal_gene_descr
        )
        if not normalized_gene_response[0] and normalized_gene_response[1]:
            return None, [normalized_gene_response[1]]

        return (
            TranscriptSegmentElement(
                transcript=genomic_data.transcript,
                exonStart=genomic_data.exon_start,
                exonStartOffset=genomic_data.exon_start_offset,
                exonEnd=genomic_data.exon_end,
                exonEndOffset=genomic_data.exon_end_offset,
                gene=normalized_gene_response[0],
                elementGenomicStart=self._sequence_location(
                    genomic_data.start,
                    genomic_data.start + 1,
                    genomic_data.chr,
                    label=genomic_data.chr,
                    seq_id_target_namespace=seq_id_target_namespace,
                )
                if genomic_data.start
                else None,
                elementGenomicEnd=self._sequence_location(
                    genomic_data.end,
                    genomic_data.end + 1,
                    genomic_data.chr,
                    label=genomic_data.chr,
                    seq_id_target_namespace=seq_id_target_namespace,
                )
                if genomic_data.end
                else None,
            ),
            None,
        )

    def gene_element(
        self, gene: str, use_minimal_gene: bool = True
    ) -> tuple[Gene | None, str | None]:
        """Create gene element

        :param str gene: Gene
        :param bool use_minimal_gene: `True` if minimal gene object
            (`id`, `gene_id`, `label`) will be used. `False` if
            gene-normalizer's gene object will be used
        :return: GeneElement, warning
        """
        return self._normalized_gene(gene, use_minimal_gene=use_minimal_gene)

    def templated_sequence_element(
        self,
        start: int,
        end: int,
        sequence_id: str,
        strand: Strand,
        label: str | None = None,
        add_location_id: bool = False,
        residue_mode: ResidueMode = ResidueMode.RESIDUE,
        seq_id_target_namespace: str | None = None,
    ) -> TemplatedSequenceElement:
        """Create templated sequence element

        :param start: Genomic start
        :param end: Genomic end
        :param sequence_id: Chromosome accession for sequence
        :param strand: Strand
        :param label: Label for genomic location
        :param add_location_id: ``True`` if ``location_id`` will be added to ``region``.
            ``False`` otherwise.
        :param residue_mode: Determines coordinate base used. Must be one of ``residue``
            or ``inter-residue``.
        :param seq_id_target_namespace: If want to use digest for ``sequence_id``, set
            this to the namespace you want the digest for. Otherwise, leave as ``None``.
        :return: Templated Sequence Element
        """
        if residue_mode == ResidueMode.RESIDUE:
            start -= 1

        region = self._sequence_location(
            start,
            end,
            sequence_id,
            label=label,
            seq_id_target_namespace=seq_id_target_namespace,
        )

        if add_location_id:
            location_id = self._location_id(region.location.model_dump())
            region.location_id = location_id

        return TemplatedSequenceElement(region=region, strand=strand)

    @staticmethod
    def linker_element(
        sequence: str,
    ) -> tuple[LinkerElement | None, str | None]:
        """Create linker element

        :param sequence: Sequence
        :return: Tuple containing a complete Linker element and None if
            successful, or a None value and warning message if unsuccessful
        """
        try:
            upper_seq = sequence.upper()
            seq = SequenceString(upper_seq)
            linker_sequence = LiteralSequenceExpression(sequence=seq)
            return LinkerElement(linkerSequence=linker_sequence), None
        except ValidationError as e:
            msg = str(e)
            _logger.warning(msg)
            return None, msg

    @staticmethod
    def multiple_possible_genes_element() -> MultiplePossibleGenesElement:
        """Create a MultiplePossibleGenesElement.

        :return: MultiplePossibleGenesElement
        """
        return MultiplePossibleGenesElement()

    @staticmethod
    def unknown_gene_element() -> UnknownGeneElement:
        """Create unknown gene element

        :return: Unknown Gene element
        """
        return UnknownGeneElement()

    def functional_domain(
        self,
        status: DomainStatus,
        name: str,
        functional_domain_id: CURIE,
        gene: str,
        sequence_id: str,
        start: int,
        end: int,
        use_minimal_gene_descr: bool = True,
        seq_id_target_namespace: str | None = None,
    ) -> tuple[FunctionalDomain | None, str | None]:
        """Build functional domain instance.

        :param status: Status for domain.  Must be either ``lost`` or ``preserved``
        :param name: Domain name
        :param functional_domain_id: Domain ID
        :param gene: Gene
        :param sequence_id: protein sequence on which provided coordinates are located
        :param start: start position on sequence
        :param end: end position on sequence
        :param use_minimal_gene_descr: ``True`` if minimal gene object (``id``,
            ``gene_id``, ``label``) will be used. ``False`` if gene-normalizer's gene
            object will be used
        :param seq_id_target_namespace: If want to use digest for ``sequence_id``, set
            this to the namespace you want the digest for. Otherwise, leave as ``None``.
        :return: Tuple with FunctionalDomain and None value for warnings if
            successful, or a None value and warning message if unsuccessful
        """
        sequence_id_lower = sequence_id.lower()
        if not (sequence_id_lower.startswith("np_")) or (
            sequence_id_lower.startswith("ensp")
        ):
            msg = "Sequence_id must be a protein accession."
            _logger.warning(msg)
            return None, msg

        seq, warning = self.cool_seq_tool.seqrepo_access.get_reference_sequence(
            sequence_id, start, end
        )

        if not seq:
            return None, warning

        gene_descr, warning = self._normalized_gene(
            gene, use_minimal_gene_descr=use_minimal_gene_descr
        )
        if not gene_descr:
            return None, warning

        loc_descr = self._sequence_location(
            start, end, sequence_id, seq_id_target_namespace=seq_id_target_namespace
        )

        try:
            return (
                FunctionalDomain(
                    id=functional_domain_id,
                    label=name,
                    status=status,
                    associatedGene=gene_descr,
                    sequence_location=loc_descr,
                ),
                None,
            )
        except ValidationError as e:
            msg = str(e)
            _logger.warning(msg)
            return None, msg

    def regulatory_element(
        self,
        regulatory_class: RegulatoryClass,
        gene: str,
        use_minimal_gene_descr: bool = True,
    ) -> tuple[RegulatoryElement | None, str | None]:
        """Create RegulatoryElement
        :param regulatory_class: one of {"promoter", "enhancer"}
        :param gene: gene term to fetch normalized gene object for
        :param use_minimal_gene_descr: whether to use the minimal gene object
        :return: Tuple with RegulatoryElement instance and None value for warnings if
            successful, or a None value and warning message if unsuccessful
        """
        gene_descr, warning = self._normalized_gene(
            gene, use_minimal_gene_descr=use_minimal_gene_descr
        )
        if not gene_descr:
            return None, warning

        try:
            return (
                RegulatoryElement(
                    regulatoryClass=regulatory_class, associatedGene=gene_descr
                ),
                None,
            )
        except ValidationError as e:
            msg = str(e)
            _logger.warning(msg)
            return None, msg

    def _sequence_location(
        self,
        start: int,
        end: int,
        sequence_id: str,
        label: str | None = None,
        seq_id_target_namespace: str | None = None,
    ) -> SequenceLocation:
        """Create sequence location

        :param start: Start position
        :param end: End position
        :param sequence_id: Accession for sequence
        :param label: label for location. If ``None``, ``sequence_id`` will be used as
            Sequence Location's ``id`` Else, label will be used as Sequence Location's
            ``id``.
        :param seq_id_target_namespace: If want to use digest for ``sequence_id``, set
            this to the namespace you want the digest for. Otherwise, leave as ``None``.
        """
        try:
            sequence_id = coerce_namespace(sequence_id)
        except ValueError:
            if not re.match(CURIE.__metadata__[0].pattern, sequence_id):
                sequence_id = f"sequence.id:{sequence_id}"

        if seq_id_target_namespace:
            try:
                seq_id = translate_identifier(
                    self.seqrepo, sequence_id, target_namespace=seq_id_target_namespace
                )
            except IDTranslationException:
                _logger.warning(
                    "Unable to translate %s using %s as the target namespace",
                    sequence_id,
                    seq_id_target_namespace,
                )
            else:
                sequence_id = seq_id

        return SequenceLocation(
            id=sequence_id,
            label=label,
            start=start,
            end=end,
        )

    def add_additional_fields(
        self,
        fusion: Fusion,
        add_all: bool = True,
        fields: list[AdditionalFields] | None = None,
        target_namespace: str = "ga4gh",
    ) -> Fusion:
        """Add additional fields to Fusion object.

        Possible fields are shown in ``AdditionalFields``

        :param fusion: A valid Fusion object
        :param add_all: ``True`` if all additional fields  will be added in fusion
            object. ``False`` if only select fields will be provided. If set to
            ``True``, will always take precedence over ``fields``.
        :param fields: Select fields that will be set. Must be a subset of
            ``AdditionalFields``
        :param target_namespace: The namespace of identifiers to return for
            ``sequence_id``. Default is ``ga4gh``
        :return: Updated fusion with specified fields set
        """
        if add_all:
            self.add_translated_sequence_id(fusion, target_namespace)
            self.add_location_id(fusion)
        else:
            if fields:
                for field in fields:
                    if field == AdditionalFields.SEQUENCE_ID.value:
                        self.add_translated_sequence_id(
                            fusion, target_namespace=target_namespace
                        )
                    elif field == AdditionalFields.LOCATION_ID.value:
                        self.add_location_id(fusion)
                    else:
                        _logger.warning("Invalid field: %s", field)
        return fusion

    def add_location_id(self, fusion: Fusion) -> Fusion:
        """Add `location_id` in fusion object.

        :param fusion: A valid Fusion object.
        :return: Updated fusion with `location_id` fields set
        """
        for structural_element in fusion.structure:
            if isinstance(structural_element, TemplatedSequenceElement):
                location = structural_element.region.location
                location_id = self._location_id(location.model_dump())
                structural_element.region.location_id = location_id
            elif isinstance(structural_element, TranscriptSegmentElement):
                for element_genomic in [
                    structural_element.elementGenomicStart,
                    structural_element.elementGenomicEnd,
                ]:
                    if element_genomic:
                        location = element_genomic.location
                        if location.type == SequenceLocation:
                            location_id = self._location_id(location.model_dump())
                            element_genomic.location_id = location_id
        if isinstance(fusion, CategoricalFusion) and fusion.criticalFunctionalDomains:
            for domain in fusion.criticalFunctionalDomains:
                location = domain.sequence_location.location
                location_id = self._location_id(location.model_dump())
                domain.sequence_location.location_id = location_id
        if fusion.regulatoryElement:
            element = fusion.regulatoryElement
            if element.feature_location:
                location = element.feature_location
                if location.type == SequenceLocation:
                    location_id = self._location_id(location.model_dump())
                    element.feature_location.location_id = location_id
        return fusion

    @staticmethod
    def _location_id(location: dict) -> CURIE:
        """Return GA4GH digest for location

        :param location: VRS Location represented as a dict
        :return: GA4GH digest
        """
        return ga4gh_identify(models.Location(**location))

    def add_translated_sequence_id(
        self, fusion: Fusion, target_namespace: str = "ga4gh"
    ) -> Fusion:
        """Translate sequence_ids in fusion object.

        :param fusion: A valid Fusion object
        :param target_namespace: ID namespace to translate sequence IDs to
        :return: Updated fusion with ``sequence_id`` fields set
        """
        for element in fusion.structure:
            if isinstance(element, TemplatedSequenceElement):
                location = element.region.location
                if location.type == SequenceLocation:
                    try:
                        new_id = translate_identifier(
                            self.seqrepo, location.sequence_id, target_namespace
                        )
                    except IDTranslationException:
                        pass
                    else:
                        element.region.location.sequence_id = new_id
            elif isinstance(element, TranscriptSegmentElement):
                for loc_descr in [
                    element.elementGenomicStart,
                    element.elementGenomicEnd,
                ]:
                    if loc_descr:
                        location = loc_descr.location
                        if location.type == SequenceLocation:
                            try:
                                new_id = translate_identifier(
                                    self.seqrepo, location.sequence_id, target_namespace
                                )
                            except IDTranslationException:
                                continue
                            loc_descr.location.sequence_id = new_id
        if fusion.type == "CategoricalFusion" and fusion.criticalFunctionalDomains:
            # TODO: unreachable code?
            for domain in fusion.criticalFunctionalDomains:
                if (
                    domain.sequence_location
                    and domain.sequence_location.location
                    and (domain.sequence_location.location.type == "SequenceLocation")
                ):
                    try:
                        new_id = translate_identifier(
                            self.seqrepo,
                            domain.sequence_location.location.sequence_id,
                            target_namespace,
                        )
                    except IDTranslationException:
                        continue
                    domain.sequence_location.location.sequence_id = new_id
        return fusion

    def _normalized_gene(
        self, query: str, use_minimal_gene: bool
    ) -> tuple[Gene | None, str | None]:
        """Return gene from normalized response.

        :param query: Gene query
        :param use_minimal_gene: bool Use minimal gene representation (id and label only)
        :return: Tuple with gene and None value for warnings if
            successful, and None value with warning string if unsuccessful
        """
        gene_norm_resp = self.gene_normalizer.normalize(query)
        if gene_norm_resp.match_type:
            gene = gene_norm_resp.gene
            if not use_minimal_gene:
                return gene, None
            return Gene(id=gene_norm_resp.normalized_id, label=gene.label), None
        return None, f"gene-normalizer unable to normalize {query}"

    def generate_nomenclature(self, fusion: Fusion) -> str:
        """Generate human-readable nomenclature describing provided fusion

        :param fusion: a valid fusion
        :return: string summarizing fusion in human-readable way per VICC fusion
            curation nomenclature
        """
        parts = []
        element_genes = []
        if fusion.regulatoryElement:
            parts.append(
                reg_element_nomenclature(fusion.regulatoryElement, self.seqrepo)
            )
        for element in fusion.structure:
            if isinstance(element, MultiplePossibleGenesElement):
                parts.append("v")
            elif isinstance(element, UnknownGeneElement):
                parts.append("?")
            elif isinstance(element, LinkerElement):
                parts.append(element.linker_sequence.sequence)
            elif isinstance(element, TranscriptSegmentElement):
                if not any(
                    [gene == element.gene.label for gene in element_genes]  # noqa: C419
                ):
                    parts.append(tx_segment_nomenclature(element))
            elif isinstance(element, TemplatedSequenceElement):
                parts.append(templated_seq_nomenclature(element, self.seqrepo))
            elif isinstance(element, GeneElement):
                if not any(
                    [gene == element.gene.label for gene in element_genes]  # noqa: C419
                ):
                    parts.append(gene_nomenclature(element))
            else:
                raise ValueError
        if (
            isinstance(fusion, AssayedFusion)
            and fusion.assay
            and fusion.assay.fusion_detection == Evidence.INFERRED
        ):
            divider = "(::)"
        else:
            divider = "::"
        return divider.join(parts)
