"""Module for testing the fusion model."""

import pytest
from cool_seq_tool.schemas import Strand
from pydantic import ValidationError

from fusor.models import (
    AbstractFusion,
    AnchoredReads,
    Assay,
    AssayedFusion,
    BreakpointCoverage,
    CategoricalFusion,
    CausativeEvent,
    ContigSequence,
    EventType,
    FunctionalDomain,
    GeneElement,
    LinkerElement,
    MultiplePossibleGenesElement,
    ReadData,
    RegulatoryElement,
    SpanningReads,
    SplitReads,
    TemplatedSequenceElement,
    TranscriptSegmentElement,
    UnknownGeneElement,
)


@pytest.fixture(scope="module")
def gene_examples():
    """Provide possible gene input."""
    return [
        {
            "primaryCoding": {
                "id": "hgnc:9339",
                "code": "HGNC:9339",
                "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
            },
            "name": "G1",
        },
        {
            "primaryCoding": {
                "id": "hgnc:76",
                "code": "HGNC:76",
                "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
            },
            "name": "ABL1",
        },
        {
            "primaryCoding": {
                "id": "hgnc:1014",
                "code": "HGNC:1014",
                "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
            },
            "name": "BCR1",
        },
        {
            "primaryCoding": {
                "id": "hgnc:8031",
                "code": "HGNC:8031",
                "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
            },
            "name": "NTRK1",
        },
        {
            "primaryCoding": {
                "id": "hgnc:1837",
                "code": "HGNC:1837",
                "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
            },
            "name": "ALK",
        },
        {
            "primaryCoding": {
                "id": "hgnc:16262",
                "code": "HGNC:16262",
                "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
            },
            "name": "YAP1",
        },
        # alternate structure
        {
            "primaryCoding": {
                "id": "hgnc:1097",
                "code": "HGNC:1097",
                "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
            },
            "type": "Gene",
            "name": "BRAF",
        },
    ]


@pytest.fixture(scope="module")
def sequence_locations():
    """Provide possible sequence_location input."""
    return [
        {
            "id": "ga4gh:SL.-xC3omZDIKZEuotbbHWQMTC8sS3nOxTb",
            "type": "SequenceLocation",
            "sequenceReference": {
                "id": "refseq:NC_000001.11",
                "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                "type": "SequenceReference",
            },
            "start": 15455,
            "end": 15456,
        },
        {
            "id": "ga4gh:SL.-xC3omZDIKZEuotbbHWQMTC8sS3nOxTb",
            "type": "SequenceLocation",
            "sequenceReference": {
                "id": "refseq:NC_000001.11",
                "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                "type": "SequenceReference",
            },
            "start": 15565,
            "end": 15566,
        },
        {
            "id": "ga4gh:SL.PPQ-aYd6dsSj7ulUEeqK8xZJP-yPrfdP",
            "type": "SequenceLocation",
            "sequenceReference": {
                "id": "refseq:NC_000012.12",
                "refgetAccession": "SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl",
                "type": "SequenceReference",
            },
            "start": 1,
            "end": 2,
        },
        {
            "id": "ga4gh:SL.OBeSv2B0pURlocL7viFiRwajew_GYGqN",
            "type": "SequenceLocation",
            "sequenceReference": {
                "id": "refseq:NC_000012.12",
                "refgetAccession": "SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl",
                "type": "SequenceReference",
            },
            "start": 2,
            "end": 3,
        },
        {
            "id": "ga4gh:SL.OBeSv2B0pURlocL7viFiRwajew_GYGqN",
            "type": "SequenceLocation",
            "sequenceReference": {
                "id": "refseq:NC_000012.12",
                "refgetAccession": "SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl",
                "type": "SequenceReference",
            },
            "start": 1,
            "end": 3,
        },
        {
            "id": "ga4gh:SL.-xC3omZDIKZEuotbbHWQMTC8sS3nOxTb",
            "type": "SequenceLocation",
            "sequenceReference": {
                "id": "refseq:NC_000001.11",
                "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                "type": "SequenceReference",
            },
            "start": 15455,
            "end": 15566,
        },
        {
            "id": "ga4gh:SL.VJLxl42yYoa-0ZMa8dfakhZfcP0nWgpl",
            "type": "SequenceLocation",
            "sequenceReference": {
                "id": "refseq:NP_001123617.1",
                "refgetAccession": "SQ.sv5egNzqN5koJQH6w0M4tIK9tEDEfJl7",
                "type": "SequenceReference",
            },
            "start": 171,
            "end": 204,
        },
        {
            "id": "ga4gh:SL.fZQW-qJwKlrVdae-idN_XXee5VTfEOgA",
            "type": "SequenceLocation",
            "sequenceReference": {
                "id": "refseq:NP_002520.2",
                "refgetAccession": "SQ.vJvm06Wl5J7DXHynR9ksW7IK3_3jlFK6",
                "type": "SequenceReference",
            },
            "start": 510,
            "end": 781,
        },
    ]


@pytest.fixture(scope="module")
def functional_domains(gene_examples, sequence_locations):
    """Provide possible functional_domains input."""
    return [
        {
            "type": "FunctionalDomain",
            "status": "preserved",
            "label": "WW domain",
            "id": "interpro:IPR001202",
            "associatedGene": gene_examples[5],
            "sequenceLocation": sequence_locations[6],
        },
        {
            "status": "lost",
            "label": "Tyrosine-protein kinase, catalytic domain",
            "id": "interpro:IPR020635",
            "associatedGene": gene_examples[3],
            "sequenceLocation": sequence_locations[7],
        },
    ]


@pytest.fixture(scope="module")
def transcript_segments(sequence_locations, gene_examples):
    """Provide possible transcript_segment input."""
    return [
        {
            "transcript": "refseq:NM_152263.3",
            "strand": -1,
            "exonStart": 1,
            "exonStartOffset": -9,
            "exonEnd": 8,
            "exonEndOffset": 7,
            "gene": gene_examples[0],
            "elementGenomicStart": sequence_locations[2],
            "elementGenomicEnd": sequence_locations[3],
            "coverage": BreakpointCoverage(fragmentCoverage=100),
            "anchoredReads": AnchoredReads(reads=85),
        },
        {
            "type": "TranscriptSegmentElement",
            "transcript": "refseq:NM_034348.3",
            "strand": 1,
            "exonStart": 1,
            "exonEnd": 8,
            "gene": gene_examples[3],
            "elementGenomicStart": sequence_locations[0],
            "elementGenomicEnd": sequence_locations[1],
        },
        {
            "type": "TranscriptSegmentElement",
            "transcript": "refseq:NM_938439.4",
            "strand": 1,
            "exonStart": 7,
            "exonEnd": 14,
            "exonEndOffset": -5,
            "gene": gene_examples[4],
            "elementGenomicStart": sequence_locations[0],
            "elementGenomicEnd": sequence_locations[1],
        },
        {
            "type": "TranscriptSegmentElement",
            "transcript": "refseq:NM_938439.4",
            "strand": 1,
            "exonStart": 7,
            "gene": gene_examples[4],
            "elementGenomicStart": sequence_locations[0],
        },
    ]


@pytest.fixture(scope="module")
def gene_elements(gene_examples):
    """Provide possible gene element input data."""
    return [
        {
            "type": "GeneElement",
            "gene": gene_examples[1],
        },
        {"type": "GeneElement", "gene": gene_examples[0]},
        {"type"},
    ]


@pytest.fixture(scope="module")
def templated_sequence_elements(sequence_locations):
    """Provide possible templated sequence element input data."""
    return [
        {
            "type": "TemplatedSequenceElement",
            "strand": 1,
            "region": sequence_locations[5],
        },
        {
            "type": "TemplatedSequenceElement",
            "strand": -1,
            "region": sequence_locations[4],
        },
    ]


@pytest.fixture(scope="module")
def literal_sequence_expressions():
    """Provide possible LiteralSequenceExpression input data"""
    return [
        {
            "id": "sequence:ACGT",
            "type": "LiteralSequenceExpression",
            "sequence": "ACGT",
        },
        {
            "id": "sequence:T",
            "type": "LiteralSequenceExpression",
            "sequence": "T",
        },
        {
            "id": "sequence:actgu",
            "type": "LiteralSequenceExpression",
            "sequence": "ACTGU",
        },
    ]


@pytest.fixture(scope="module")
def linkers(literal_sequence_expressions):
    """Provide possible linker element input data."""
    return [
        {
            "type": "LinkerSequenceElement",
            "linkerSequence": literal_sequence_expressions[0],
        },
        {
            "type": "LinkerSequenceElement",
            "linkerSequence": literal_sequence_expressions[1],
        },
        {
            "type": "LinkerSequenceElement",
            "linkerSequence": literal_sequence_expressions[2],
        },
    ]


@pytest.fixture(scope="module")
def unknown_element():
    """Provide UnknownGene element."""
    return {"type": "UnknownGeneElement"}


@pytest.fixture(scope="module")
def regulatory_elements(gene_examples):
    """Provide possible regulatory_element input data."""
    return [{"regulatoryClass": "promoter", "associatedGene": gene_examples[0]}]


def check_validation_error(exc_info, expected_msg: str, index: int = 0):
    """Check ValidationError instance for expected message.

    :param ExceptionInfo exc_info: ValidationError instance raised and captured
    by pytest.
    :param str expected_msg: message expected to be provided by error
    :param int index: optional index (if multiple errors are raised)
    :return: None, but may raise AssertionError if incorrect behavior found.
    """
    assert exc_info.value.errors()[index]["msg"] == expected_msg


def test_functional_domain(functional_domains, gene_examples):
    """Test FunctionalDomain object initializes correctly"""
    test_domain = FunctionalDomain(**functional_domains[0])
    assert test_domain.type == "FunctionalDomain"
    assert test_domain.status == "preserved"
    assert test_domain.label == "WW domain"
    assert test_domain.id == "interpro:IPR001202"
    assert test_domain.associatedGene.primaryCoding.id == "hgnc:16262"
    assert test_domain.associatedGene.name == "YAP1"
    test_loc = test_domain.sequenceLocation
    assert "ga4gh:SL" in test_loc.id
    assert test_loc.type == "SequenceLocation"
    assert test_loc.start == 171
    assert test_loc.end == 204
    test_ref = test_loc.sequenceReference
    assert test_ref.id == "refseq:NP_001123617.1"
    assert "SQ." in test_ref.refgetAccession
    assert test_ref.type == "SequenceReference"

    test_domain = FunctionalDomain(**functional_domains[1])
    assert test_domain.type == "FunctionalDomain"
    assert test_domain.status == "lost"
    assert test_domain.label == "Tyrosine-protein kinase, catalytic domain"
    assert test_domain.id == "interpro:IPR020635"
    assert test_domain.associatedGene.primaryCoding.id == "hgnc:8031"
    assert test_domain.associatedGene.name == "NTRK1"
    test_loc = test_domain.sequenceLocation
    assert "ga4gh:SL" in test_loc.id
    assert test_loc.type == "SequenceLocation"
    assert test_loc.start == 510
    assert test_loc.end == 781
    test_ref = test_loc.sequenceReference
    assert test_ref.id == "refseq:NP_002520.2"
    assert "SQ." in test_ref.refgetAccession
    assert test_ref.type == "SequenceReference"

    # test status string
    with pytest.raises(ValidationError) as exc_info:
        FunctionalDomain(
            status="gained",
            name="tyrosine kinase catalytic domain",
            id="interpro:IPR020635",
            associatedGene=gene_examples[0],
        )
    msg = "Input should be 'lost' or 'preserved'"
    check_validation_error(exc_info, msg)

    # test domain ID CURIE requirement
    with pytest.raises(ValidationError) as exc_info:
        FunctionalDomain(
            status="lost",
            label="tyrosine kinase catalytic domain",
            id="interpro_IPR020635",
            associatedGene=gene_examples[0],
        )
    msg = "String should match pattern '^\\w[^:]*:.+$'"
    check_validation_error(exc_info, msg)


def test_transcript_segment_element(transcript_segments):
    """Test TranscriptSegmentElement object initializes correctly"""
    test_element = TranscriptSegmentElement(**transcript_segments[0])
    assert test_element.transcript == "refseq:NM_152263.3"
    assert test_element.strand == -1
    assert test_element.exonStart == 1
    assert test_element.exonStartOffset == -9
    assert test_element.exonEnd == 8
    assert test_element.exonEndOffset == 7
    assert test_element.gene.primaryCoding.id == "hgnc:9339"
    assert test_element.gene.name == "G1"
    test_region_start = test_element.elementGenomicStart
    assert test_region_start.type == "SequenceLocation"
    test_region_end = test_element.elementGenomicEnd
    assert test_region_end.type == "SequenceLocation"
    assert test_element.coverage.fragmentCoverage == 100
    assert test_element.anchoredReads.reads == 85

    test_element = TranscriptSegmentElement(**transcript_segments[3])
    assert test_element.transcript == "refseq:NM_938439.4"
    assert test_element.strand == 1
    assert test_element.exonStart == 7
    assert test_element.exonStartOffset == 0
    assert test_element.exonEnd is None
    assert test_element.exonEndOffset is None
    assert test_element.coverage is None
    assert test_element.anchoredReads is None

    # check CURIE requirement
    with pytest.raises(ValidationError) as exc_info:
        TranscriptSegmentElement(
            transcript="NM_152263.3",
            strand="-1",
            exonStart="1",
            exonStartOffset="-9",
            exonEnd="8",
            exonEndOffset="7",
            gene={
                "primaryCoding": {"id": "test:1", "code": "test:1", "system": "test"},
                "name": "G1",
            },
            elementGenomicStart={
                "location": {
                    "species_id": "taxonomy:9606",
                    "chr": "12",
                    "interval": {"start": "p12.1", "end": "p12.1"},
                }
            },
            elementGenomicEnd={
                "location": {
                    "species_id": "taxonomy:9606",
                    "chr": "12",
                    "interval": {"start": "p12.2", "end": "p12.2"},
                }
            },
        )
    msg = "String should match pattern '^\\w[^:]*:.+$'"
    check_validation_error(exc_info, msg)

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert TranscriptSegmentElement(
            type="TemplatedSequenceElement",
            transcript="NM_152263.3",
            strand="-1",
            exonStart="1",
            exonStartOffset="-9",
            exonEnd="8",
            exonEndOffset="7",
            gene={
                "primaryCode": "test:1",
                "name": "G1",
            },
            elementGenomicStart={
                "location": {
                    "species_id": "taxonomy:9606",
                    "chr": "12",
                    "interval": {"start": "p12.1", "end": "p12.2"},
                }
            },
            elementGenomicEnd={
                "location": {
                    "species_id": "taxonomy:9606",
                    "chr": "12",
                    "interval": {"start": "p12.2", "end": "p12.2"},
                }
            },
        )
    msg = "Input should be <FUSORTypes.TRANSCRIPT_SEGMENT_ELEMENT: 'TranscriptSegmentElement'>"
    check_validation_error(exc_info, msg)

    # test element required
    with pytest.raises(ValidationError) as exc_info:
        assert TranscriptSegmentElement(
            element_type="templated_sequence",
            transcript="NM_152263.3",
            strand="-1",
            exonStart="1",
            exonStartOffset="-9",
            gene={
                "primaryCoding": {"id": "test:1", "code": "test:1", "system": "test"},
                "name": "G1",
            },
        )
    msg = "Value error, Must give `elementGenomicStart` if `exonStart` is given"
    check_validation_error(exc_info, msg)

    # Neither exonStart or exonEnd given
    with pytest.raises(ValidationError) as exc_info:
        assert TranscriptSegmentElement(
            type="TranscriptSegmentElement",
            transcript="NM_152263.3",
            strand="-1",
            exonStartOffset="-9",
            exonEndOffset="7",
            gene={
                "primaryCoding": {"id": "test:1", "code": "test:1", "system": "test"},
                "name": "G1",
            },
            elementGenomicStart={
                "location": {
                    "species_id": "taxonomy:9606",
                    "chr": "12",
                    "interval": {"start": "p12.1", "end": "p12.2"},
                }
            },
            elementGenomicEnd={
                "location": {
                    "species_id": "taxonomy:9606",
                    "chr": "12",
                    "interval": {"start": "p12.2", "end": "p12.2"},
                }
            },
        )
    msg = "Value error, Must give values for either `exonStart`, `exonEnd`, or both"
    check_validation_error(exc_info, msg)


def test_linker_element(linkers):
    """Test Linker object initializes correctly"""

    def check_linker(actual, expected_id, expected_sequence):
        assert actual.type == "LinkerSequenceElement"
        assert actual.linkerSequence.id == expected_id
        assert actual.linkerSequence.sequence.root == expected_sequence
        assert actual.linkerSequence.type == "LiteralSequenceExpression"

    for args in (
        (LinkerElement(**linkers[0]), "sequence:ACGT", "ACGT"),
        (LinkerElement(**linkers[1]), "sequence:T", "T"),
        (LinkerElement(**linkers[2]), "sequence:actgu", "ACTGU"),
    ):
        check_linker(*args)

    # check base validation
    with pytest.raises(ValidationError) as exc_info:
        LinkerElement(linkerSequence={"id": "sequence:ACT1", "sequence": "ACT1"})
    msg = "String should match pattern '^[A-Z*\\-]*$'"
    check_validation_error(exc_info, msg)

    # check valid literal sequence expression
    with pytest.raises(ValidationError) as exc_info:
        LinkerElement(linkerSequence={"id": "sequence:actgu", "sequence": "actgu"})
    msg = "String should match pattern '^[A-Z*\\-]*$'"
    check_validation_error(exc_info, msg)

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert LinkerElement(
            type="TemplatedSequenceElement",
            linkerSequence={"id": "sequence:ATG", "sequence": "ATG"},
        )
    msg = (
        "Input should be <FUSORTypes.LINKER_SEQUENCE_ELEMENT: 'LinkerSequenceElement'>"
    )
    check_validation_error(exc_info, msg)

    # test no extras
    with pytest.raises(ValidationError) as exc_info:
        assert LinkerElement(
            type="LinkerSequenceElement",
            linkerSequence={"id": "sequence:G", "sequence": "G"},
            bonus_value="bonus",
        )
    msg = "Extra inputs are not permitted"
    check_validation_error(exc_info, msg)


def test_genomic_region_element(templated_sequence_elements, sequence_locations):
    """Test that TemplatedSequenceElement initializes correctly."""

    def assert_genomic_region_test_element(test):
        """Assert that test templated_sequence_elements[0] data matches
        expected values.
        """
        assert test.type == "TemplatedSequenceElement"
        assert test.strand == Strand.POSITIVE
        assert "ga4gh:SL" in test.region.id
        assert test.region.type == "SequenceLocation"
        test_ref = test.region.sequenceReference
        assert "refseq:" in test_ref.id
        assert "SQ." in test_ref.refgetAccession

    test_element = TemplatedSequenceElement(**templated_sequence_elements[0])
    assert_genomic_region_test_element(test_element)

    with pytest.raises(ValidationError) as exc_info:
        TemplatedSequenceElement(
            region={
                "start": 39408,
                "end": 39414,
                "sequenceReference": {
                    "id": "refseq:NC_000012.12",
                    "refgetAccession": "SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl",
                },
            },
        )
    msg = "Field required"
    check_validation_error(exc_info, msg)

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert TemplatedSequenceElement(
            type="GeneElement", region=sequence_locations[0], strand=Strand.POSITIVE
        )
    msg = "Input should be <FUSORTypes.TEMPLATED_SEQUENCE_ELEMENT: 'TemplatedSequenceElement'>"
    check_validation_error(exc_info, msg)


def test_gene_element(gene_examples):
    """Test that Gene Element initializes correctly."""
    test_element = GeneElement(gene=gene_examples[0])
    assert test_element.type == "GeneElement"
    assert test_element.gene.primaryCoding.id == "hgnc:9339"
    assert test_element.gene.name == "G1"

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert GeneElement(type="UnknownGeneElement", gene=gene_examples[0])
    msg = "Input should be <FUSORTypes.GENE_ELEMENT: 'GeneElement'>"
    check_validation_error(exc_info, msg)


def test_unknown_gene_element():
    """Test that unknown_gene element initializes correctly."""
    test_element = UnknownGeneElement()
    assert test_element.type == "UnknownGeneElement"

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert UnknownGeneElement(type="gene")
    msg = "Input should be <FUSORTypes.UNKNOWN_GENE_ELEMENT: 'UnknownGeneElement'>"
    check_validation_error(exc_info, msg)


def test_mult_gene_element():
    """Test that mult_gene_element initializes correctly."""
    test_element = MultiplePossibleGenesElement()
    assert test_element.type == "MultiplePossibleGenesElement"

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert MultiplePossibleGenesElement(type="unknown_gene")
    msg = "Input should be <FUSORTypes.MULTIPLE_POSSIBLE_GENES_ELEMENT: 'MultiplePossibleGenesElement'>"
    check_validation_error(exc_info, msg)


def test_coverage():
    """Test that BreakpointCoverage class initializes correctly"""
    test_coverage = BreakpointCoverage(fragmentCoverage=100)
    assert test_coverage.fragmentCoverage == 100

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert BreakpointCoverage(type="coverage")
    msg = "Input should be <FUSORTypes.BREAKPOINT_COVERAGE: 'BreakpointCoverage'>"
    check_validation_error(exc_info, msg)


def test_contig():
    """Test that Contig class initializes correctly"""
    test_contig = ContigSequence(contig="GTATACTATGATCAGT")
    assert test_contig.contig == "GTATACTATGATCAGT"

    test_contig = ContigSequence(contig="GTATACTATGATCAGT|ATGATCATGAT")
    assert test_contig.contig == "GTATACTATGATCAGT|ATGATCATGAT"

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert ContigSequence(type="contig")
    msg = "Input should be <FUSORTypes.CONTIG_SEQUENCE: 'ContigSequence'>"
    check_validation_error(exc_info, msg)

    # test invalid input
    with pytest.raises(ValidationError) as exc_info:
        ContigSequence(contig="1212341|ATGATCATGAT")
    msg = "String should match pattern '^(?:[^A-Za-z0-9]|[ACTGactg])*$'"
    check_validation_error(exc_info, msg)


def test_anchored_reads():
    """Test that AnchoredReads class initializes correctly"""
    test_anchored_reads = AnchoredReads(reads=100)
    assert test_anchored_reads.reads == 100

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert AnchoredReads(type="anchoredreads")
    msg = "Input should be <FUSORTypes.ANCHORED_READS: 'AnchoredReads'>"
    check_validation_error(exc_info, msg)


def test_split_reads():
    """Test that SplitReads class initializes correctly"""
    test_split_reads = SplitReads(splitReads=97)
    assert test_split_reads.splitReads == 97

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert SplitReads(type="splitreads")
    msg = "Input should be <FUSORTypes.SPLIT_READS: 'SplitReads'>"
    check_validation_error(exc_info, msg)


def test_spanning_reads():
    """Test that SpanningReads class initializes correctly"""
    test_spanning_reads = SpanningReads(spanningReads=97)
    assert test_spanning_reads.spanningReads == 97

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert SpanningReads(type="spanningreads")
    msg = "Input should be <FUSORTypes.SPANNING_READS: 'SpanningReads'>"
    check_validation_error(exc_info, msg)


def test_read_data():
    """Test that ReadData class initializes correctly"""
    test_read_data = ReadData(
        split=SplitReads(splitReads=100), spanning=SpanningReads(spanningReads=90)
    )
    assert test_read_data.split.splitReads == 100
    assert test_read_data.spanning.spanningReads == 90

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert ReadData(type="readata")
    msg = "Input should be <FUSORTypes.READ_DATA: 'ReadData'>"
    check_validation_error(exc_info, msg)


def test_event():
    """Test Event object initializes correctly"""
    rearrangement = EventType.REARRANGEMENT
    test_event = CausativeEvent(eventType=rearrangement, eventDescription=None)
    assert test_event.eventType == rearrangement

    with pytest.raises(ValueError):  # noqa: PT011
        CausativeEvent(eventType="combination")


def test_regulatory_element(regulatory_elements, gene_examples):
    """Test RegulatoryElement object initializes correctly"""
    test_reg_elmt = RegulatoryElement(**regulatory_elements[0])
    assert test_reg_elmt.regulatoryClass.value == "promoter"
    assert test_reg_elmt.associatedGene.primaryCoding.id == "hgnc:9339"
    assert test_reg_elmt.associatedGene.name == "G1"

    # check type constraint
    with pytest.raises(ValidationError) as exc_info:
        RegulatoryElement(
            regulatoryClass="notpromoter", associatedGene=gene_examples[0]
        )
    assert exc_info.value.errors()[0]["msg"].startswith("Input should be")

    # require minimum input
    with pytest.raises(ValidationError) as exc_info:
        RegulatoryElement(regulatoryClass="enhancer")
    assert (
        exc_info.value.errors()[0]["msg"]
        == "Value error, Must set 1 of {`featureId`, `associatedGene`} and/or `featureLocation`"
    )


def test_fusion(
    functional_domains,
    transcript_segments,
    templated_sequence_elements,
    linkers,
    gene_elements,
    regulatory_elements,
    unknown_element,
):
    """Test that Fusion object initializes correctly"""
    # test valid object
    fusion = CategoricalFusion(
        readingFramePreserved=True,
        criticalFunctionalDomains=[functional_domains[0]],
        structure=[transcript_segments[1], transcript_segments[2]],
        regulatoryElement=regulatory_elements[0],
    )

    assert fusion.structure[0].transcript == "refseq:NM_034348.3"

    # check correct parsing of nested items
    fusion = CategoricalFusion(
        structure=[
            {
                "type": "GeneElement",
                "gene": {
                    "conceptType": "Gene",
                    "primaryCoding": {
                        "id": "hgnc:8031",
                        "code": "HGNC:8031",
                        "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                    },
                    "name": "NTRK1",
                },
            },
            {
                "type": "GeneElement",
                "gene": {
                    "conceptType": "Gene",
                    "primaryCoding": {
                        "id": "hgnc:76",
                        "code": "HGNC:76",
                        "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                    },
                    "name": "ABL1",
                },
            },
        ],
        regulatory_element=None,
    )
    assert fusion.structure[0].type == "GeneElement"
    assert fusion.structure[0].gene.name == "NTRK1"
    assert fusion.structure[0].gene.primaryCoding.id == "hgnc:8031"
    assert fusion.structure[1].type == "GeneElement"

    # test that non-element properties are optional
    assert CategoricalFusion(structure=[transcript_segments[1], transcript_segments[2]])

    # test variety of element types
    causative_event = {
        "type": "CausativeEvent",
        "eventType": "rearrangement",
        "eventDescription": "chr2:g.pter_8,247,756::chr11:g.15,825,273_cen_qter (der11) and chr11:g.pter_15,825,272::chr2:g.8,247,757_cen_qter (der2)",
    }
    assay = {
        "type": "Assay",
        "methodUri": "pmid:33576979",
        "assayId": "obi:OBI_0003094",
        "assayName": "fluorescence in-situ hybridization assay",
        "fusionDetection": "inferred",
    }
    assert AssayedFusion(
        type="AssayedFusion",
        structure=[
            unknown_element,
            gene_elements[0],
            transcript_segments[2],
            templated_sequence_elements[1],
            linkers[0],
        ],
        causativeEvent=causative_event,
        assay=assay,
    )
    with pytest.raises(ValidationError) as exc_info:
        assert CategoricalFusion(
            type="CategoricalFusion",
            structure=[
                {
                    "type": "LinkerSequenceElement",
                    "linkerSequence": {
                        "id": "a:b",
                        "type": "LiteralSequenceExpression",
                        "sequence": "AC",
                    },
                },
                {
                    "type": "LinkerSequenceElement",
                    "linkerSequence": {
                        "id": "a:b",
                        "type": "LiteralSequenceExpression",
                        "sequence": "AC",
                    },
                },
            ],
        )
    msg = "Value error, First structural element cannot be LinkerSequence"
    check_validation_error(exc_info, msg)

    with pytest.raises(ValidationError) as exc_info:
        assert AssayedFusion(
            type="AssayedFusion",
            structure=[
                transcript_segments[3],
                transcript_segments[1],
            ],
            causativeEvent=causative_event,
            assay=assay,
        )
    msg = "Value error, 5' TranscriptSegmentElement fusion partner must contain ending exon position"
    check_validation_error(exc_info, msg)


def test_fusion_element_count(
    functional_domains,
    regulatory_elements,
    unknown_element,
    gene_elements,
    transcript_segments,
    gene_examples,
):
    """Test fusion element count requirements."""
    # elements are mandatory
    with pytest.raises(ValidationError) as exc_info:
        assert AssayedFusion(
            functionalDomains=[functional_domains[1]],
            causativeEvent="rearrangement",
            regulatoryElement=[regulatory_elements[0]],
        )
    element_ct_msg = (
        "Value error, Fusions must contain >= 2 structural elements, or >=1 structural element "
        "and a regulatory element"
    )
    check_validation_error(exc_info, element_ct_msg)

    # must have >= 2 elements + regulatory elements
    with pytest.raises(ValidationError) as exc_info:
        assert AssayedFusion(
            structure=[unknown_element],
            causativeEvent={
                "type": "CausativeEvent",
                "eventType": "rearrangement",
                "event_description": "chr2:g.pter_8,247,756::chr11:g.15,825,273_cen_qter (der11) and chr11:g.pter_15,825,272::chr2:g.8,247,757_cen_qter (der2)",
            },
            assay={
                "type": "Assay",
                "methodUri": "pmid:33576979",
                "assayId": "obi:OBI_0003094",
                "assayName": "fluorescence in-situ hybridization assay",
                "fusionDetection": "inferred",
            },
        )
    check_validation_error(exc_info, element_ct_msg)

    # unique gene requirements
    uq_gene_error_msg = "Value error, Fusions must form a chimeric transcript from two or more genes, or a novel interaction between a rearranged regulatory element with the expressed product of a partner gene."
    with pytest.raises(ValidationError) as exc_info:
        assert CategoricalFusion(structure=[gene_elements[0], gene_elements[0]])
    check_validation_error(exc_info, uq_gene_error_msg)

    with pytest.raises(ValidationError) as exc_info:
        assert CategoricalFusion(structure=[gene_elements[1], transcript_segments[0]])
    check_validation_error(exc_info, uq_gene_error_msg)

    with pytest.raises(ValidationError) as exc_info:
        assert CategoricalFusion(
            regulatoryElement=regulatory_elements[0],
            structure=[transcript_segments[0]],
        )
    check_validation_error(exc_info, uq_gene_error_msg)

    # use alternate gene structure
    with pytest.raises(ValidationError) as exc_info:
        assert AssayedFusion(
            type="AssayedFusion",
            structure=[
                {"type": "GeneElement", "gene": gene_examples[6]},
                {"type": "GeneElement", "gene": gene_examples[6]},
            ],
            causativeEvent={
                "type": "CausativeEvent",
                "eventType": "read-through",
            },
            assay={
                "type": "Assay",
                "methodUri": "pmid:33576979",
                "assayId": "obi:OBI_0003094",
                "assayName": "fluorescence in-situ hybridization assay",
                "fusionDetection": "inferred",
            },
        )
    with pytest.raises(ValidationError) as exc_info:
        assert AssayedFusion(
            type="AssayedFusion",
            structure=[
                {"type": "GeneElement", "gene": gene_examples[6]},
            ],
            regulatoryElement={
                "type": "RegulatoryElement",
                "regulatory_class": "enhancer",
                "feature_id": "EH111111111",
                "associatedGene": gene_examples[6],
            },
            causativeEvent={
                "type": "CausativeEvent",
                "eventType": "read-through",
            },
            assay={
                "type": "Assay",
                "methodUri": "pmid:33576979",
                "assayId": "obi:OBI_0003094",
                "assayName": "fluorescence in-situ hybridization assay",
                "fusionDetection": "inferred",
            },
        )


def test_fusion_abstraction_validator(transcript_segments, linkers):
    """Test that instantiation of abstract fusion fails."""
    # can't create base fusion
    with pytest.raises(ValidationError) as exc_info:
        assert AbstractFusion(structure=[transcript_segments[2], linkers[0]])
    check_validation_error(
        exc_info, "Value error, Cannot instantiate Fusion abstract class"
    )


def test_file_examples():
    """Test example JSON files."""
    # if this loads, then Pydantic validation was successful
    import fusor.examples as _  # noqa: F401


def test_model_examples():
    """Test example objects as provided in Pydantic config classes"""
    models = [
        FunctionalDomain,
        TranscriptSegmentElement,
        LinkerElement,
        TemplatedSequenceElement,
        GeneElement,
        UnknownGeneElement,
        MultiplePossibleGenesElement,
        RegulatoryElement,
        Assay,
        CausativeEvent,
        AssayedFusion,
        CategoricalFusion,
    ]
    for model in models:
        schema = model.model_config["json_schema_extra"]
        if "example" in schema:
            model(**schema["example"])
