"""Module for testing the FUSOR class."""

import copy

import pytest
from cool_seq_tool.schemas import Strand
from ga4gh.core.models import Coding, MappableConcept
from ga4gh.vrs.models import SequenceLocation

from fusor.exceptions import FUSORParametersException, IDTranslationException
from fusor.models import (
    AssayedFusion,
    CategoricalFusion,
    FunctionalDomain,
    GeneElement,
    LinkerElement,
    MultiplePossibleGenesElement,
    RegulatoryClass,
    RegulatoryElement,
    TemplatedSequenceElement,
    TranscriptSegmentElement,
    UnknownGeneElement,
)


@pytest.fixture(scope="module")
def braf_gene_obj_min():
    """Create minimal gene object for BRAF"""
    return MappableConcept(
        primaryCoding=Coding(
            id="hgnc:1097",
            code="HGNC:1097",
            system="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
        ),
        name="BRAF",
        conceptType="Gene",
    )


@pytest.fixture(scope="module")
def braf_gene_obj(braf_gene):
    """Create gene object for braf"""
    return MappableConcept(**braf_gene)


@pytest.fixture(scope="module")
def linker_element():
    """Create linker element test fixture."""
    params = {
        "linkerSequence": {
            "id": "fusor.sequence:act",
            "sequence": "ACT",
            "type": "LiteralSequenceExpression",
        },
        "type": "LinkerSequenceElement",
    }
    return LinkerElement(**params)


@pytest.fixture(scope="module")
def sequence_location_braf_domain():
    """Create sequence location fixture for BRAF catalytic domain"""
    params = {
        "id": "ga4gh:SL.Lm-hzZHlA8FU_cYaOtAIbMLdf4Kk-SF8",
        "type": "SequenceLocation",
        "sequenceReference": {
            "id": "refseq:NP_004324.2",
            "refgetAccession": "SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
            "type": "SequenceReference",
        },
        "start": 458,
        "end": 712,
    }
    return SequenceLocation(**params)


@pytest.fixture(scope="module")
def sequence_location_braf_ref_id_ga4gh():
    """Create sequence location fixture for BRAF catalytic domain"""
    params = {
        "id": "ga4gh:SL.Lm-hzZHlA8FU_cYaOtAIbMLdf4Kk-SF8",
        "type": "SequenceLocation",
        "sequenceReference": {
            "id": "ga4gh:SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
            "refgetAccession": "SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
            "type": "SequenceReference",
        },
        "start": 458,
        "end": 712,
    }
    return SequenceLocation(**params)


@pytest.fixture(scope="module")
def functional_domain_min(braf_gene_obj_min, sequence_location_braf_domain):
    """Create functional domain test fixture."""
    params = {
        "status": "preserved",
        "label": "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "id": "interpro:IPR001245",
        "associatedGene": braf_gene_obj_min,
        "sequenceLocation": sequence_location_braf_domain,
    }
    return FunctionalDomain(**params)


@pytest.fixture(scope="module")
def functional_domain(braf_gene_obj, sequence_location_braf_domain):
    """Create functional domain test fixture."""
    params = {
        "status": "preserved",
        "label": "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "id": "interpro:IPR001245",
        "associatedGene": braf_gene_obj,
        "sequenceLocation": sequence_location_braf_domain,
    }
    return FunctionalDomain(**params)


@pytest.fixture(scope="module")
def functional_domain_seq_id(braf_gene_obj_min, sequence_location_braf_ref_id_ga4gh):
    """Create functional domain test fixture."""
    params = {
        "status": "preserved",
        "label": "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "id": "interpro:IPR001245",
        "associatedGene": braf_gene_obj_min,
        "sequenceLocation": sequence_location_braf_ref_id_ga4gh,
    }
    return FunctionalDomain(**params)


@pytest.fixture(scope="module")
def regulatory_element(braf_gene_obj):
    """Create regulatory element test fixture."""
    params = {
        "type": "RegulatoryElement",
        "regulatoryClass": "promoter",
        "associatedGene": braf_gene_obj,
    }
    return RegulatoryElement(**params)


@pytest.fixture(scope="module")
def regulatory_element_min(braf_gene_obj_min):
    """Create regulatory element test fixture with minimal gene object."""
    params = {"regulatoryClass": "promoter", "associatedGene": braf_gene_obj_min}
    return RegulatoryElement(**params)


@pytest.fixture(scope="module")
def sequence_location_tpm3():
    """Create sequence location for TPM3 test fixture."""
    params = {
        "id": "ga4gh:SL.0cMJgKuY32ate6k95oLua6vv8JAJ4PzO",
        "type": "SequenceLocation",
        "sequenceReference": {
            "id": "NM_152263.3",
            "refgetAccession": "SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT",
            "type": "SequenceReference",
        },
        "start": 154170398,
        "end": 154170399,
    }
    return SequenceLocation(**params)


@pytest.fixture(scope="module")
def templated_sequence_element():
    """Create test fixture for templated sequence element"""
    params = {
        "type": "TemplatedSequenceElement",
        "region": {
            "id": "ga4gh:SL.U7-HtnKxK9kKI1ZINiDM_m4I6O-p4Dc9",
            "digest": "U7-HtnKxK9kKI1ZINiDM_m4I6O-p4Dc9",
            "type": "SequenceLocation",
            "sequenceReference": {
                "id": "refseq:NC_000001.11",
                "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                "type": "SequenceReference",
            },
            "start": 99,
            "end": 150,
        },
        "strand": 1,
    }
    return TemplatedSequenceElement(**params)


@pytest.fixture(scope="module")
def transcript_segment_element():
    """Create transcript segment element test fixture"""
    params = {
        "type": "TranscriptSegmentElement",
        "strand": -1,
        "exonEnd": 8,
        "exonEndOffset": 0,
        "exonStart": 1,
        "exonStartOffset": 0,
        "gene": {
            "primaryCoding": {
                "id": "hgnc:12012",
                "code": "HGNC:12012",
                "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
            },
            "name": "TPM3",
            "conceptType": "Gene",
        },
        "transcript": "refseq:NM_152263.3",
        "elementGenomicStart": {
            "id": "ga4gh:SL.Q8vkGp7_xR9vI0PQ7g1IvUUeQ4JlJG8l",
            "digest": "Q8vkGp7_xR9vI0PQ7g1IvUUeQ4JlJG8l",
            "type": "SequenceLocation",
            "sequenceReference": {
                "id": "refseq:NC_000001.11",
                "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                "type": "SequenceReference",
            },
            "end": 154192135,
        },
        "elementGenomicEnd": {
            "id": "ga4gh:SL.Lnne0bSsgjzmNkKsNnXg98FeJSrDJuLb",
            "digest": "Lnne0bSsgjzmNkKsNnXg98FeJSrDJuLb",
            "type": "SequenceLocation",
            "sequenceReference": {
                "id": "refseq:NC_000001.11",
                "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                "type": "SequenceReference",
            },
            "start": 154170399,
        },
    }
    return TranscriptSegmentElement(**params)


@pytest.fixture(scope="module")
def mane_transcript_segment_element():
    """Create transcript segment element test fixture"""
    params = {
        "type": "TranscriptSegmentElement",
        "strand": 1,
        "exonEnd": None,
        "exonEndOffset": None,
        "exonStart": 2,
        "exonStartOffset": 0,
        "gene": {
            "primaryCoding": {
                "id": "hgnc:12761",
                "code": "HGNC:12761",
                "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
            },
            "name": "WEE1",
            "conceptType": "Gene",
        },
        "transcript": "refseq:NM_003390.4",
        "elementGenomicEnd": None,
        "elementGenomicStart": {
            "id": "ga4gh:SL.Dm_Rri77OtV3-FmEmGXBjWZ2PhEzdhFT",
            "digest": "Dm_Rri77OtV3-FmEmGXBjWZ2PhEzdhFT",
            "type": "SequenceLocation",
            "sequenceReference": {
                "id": "refseq:NC_000011.10",
                "refgetAccession": "SQ.2NkFm8HK88MqeNkCgj78KidCAXgnsfV1",
                "type": "SequenceReference",
            },
            "start": 9575887,
        },
    }
    return TranscriptSegmentElement(**params)


@pytest.fixture()
def fusion_ensg_sequence_id(templated_sequence_element_ensg):
    """Create fixture using Ensemble gene ID."""
    params = {
        "type": "CategoricalFusion",
        "structure": [
            templated_sequence_element_ensg,
            {"type": "MultiplePossibleGenesElement"},
        ],
        "readingFramePreserved": True,
        "regulatoryElement": None,
    }
    return CategoricalFusion(**params)


def compare_gene_obj(actual: dict, expected: dict):
    """Test that actual and expected gene objects match."""
    assert actual["primaryCoding"] == expected["primaryCoding"]
    assert actual["name"] == expected["name"]
    assert actual["conceptType"] == expected["conceptType"]
    if expected.get("xrefs"):
        assert set(actual.get("xrefs")) == set(expected["xrefs"]), "xrefs"
    else:
        assert actual.get("xrefs") == expected.get("xrefs")
    assert "extensions" in actual
    if expected["extensions"]:
        assert len(actual["extensions"]) == len(
            expected["extensions"]
        ), "len of extensions"
        n_ext_correct = 0
        for expected_ext in expected["extensions"]:
            for actual_ext in actual["extensions"]:
                if actual_ext["name"] == expected_ext["name"]:
                    assert isinstance(actual_ext["value"], type(expected_ext["value"]))
                    if isinstance(expected_ext["value"], list) and not isinstance(
                        expected_ext["value"][0], dict
                    ):
                        assert set(actual_ext["value"]) == set(
                            expected_ext["value"]
                        ), f"{expected_ext['value']} value"
                    else:
                        assert actual_ext["value"] == expected_ext["value"]
                    assert actual_ext.get("type") == expected_ext.get("type")
                    n_ext_correct += 1
        assert n_ext_correct == len(
            expected["extensions"]
        ), "number of correct extensions"


def test__normalized_gene(fusor_instance):
    """Test that _normalized_gene works correctly."""
    # TODO: test actual response
    resp = fusor_instance._normalized_gene("BRAF")
    assert resp[0]
    assert resp[1] is None
    assert isinstance(resp[0], MappableConcept)

    resp = fusor_instance._normalized_gene("B R A F")
    assert resp[0] is None
    assert resp[1] == "gene-normalizer unable to normalize B R A F"


def test_fusion(
    fusor_instance,
    linker_element,
    templated_sequence_element,
    transcript_segment_element,
    functional_domain,
):
    """Test that fusion methods work correctly."""
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

    # infer type from properties
    f = fusor_instance.fusion(
        structure=[
            templated_sequence_element,
            linker_element,
            UnknownGeneElement(),
        ],
        causative_event=causative_event,
        assay=assay,
    )
    assert isinstance(f, AssayedFusion)
    f = fusor_instance.fusion(
        structure=[
            transcript_segment_element,
            MultiplePossibleGenesElement(),
        ],
        critical_functional_domains=[functional_domain],
    )
    assert isinstance(f, CategoricalFusion)

    # catch conflicting property args
    with pytest.raises(FUSORParametersException) as excinfo:
        f = fusor_instance.fusion(
            structure=[
                transcript_segment_element,
                UnknownGeneElement(),
            ],
            causative_event="rearrangement",
            critical_functional_domains=[functional_domain],
        )
    assert str(excinfo.value) == "Received conflicting attributes"

    # handle indeterminate type
    with pytest.raises(FUSORParametersException) as excinfo:
        f = fusor_instance.fusion(
            structure=[
                transcript_segment_element,
                templated_sequence_element,
            ]
        )
    assert str(excinfo.value) == "Unable to determine fusion type"

    # handle both type parameter options
    f = fusor_instance.fusion(
        fusion_type="AssayedFusion",
        structure=[
            templated_sequence_element,
            linker_element,
            UnknownGeneElement(),
        ],
        causative_event={
            "type": "CausativeEvent",
            "eventType": "rearrangement",
        },
        assay={
            "type": "Assay",
            "methodUri": "pmid:33576979",
            "assayId": "obi:OBI_0003094",
            "assayName": "fluorescence in-situ hybridization assay",
            "fusionDetection": "inferred",
        },
    )
    assert isinstance(f, AssayedFusion)
    f = fusor_instance.fusion(
        type="CategoricalFusion",
        structure=[
            transcript_segment_element,
            MultiplePossibleGenesElement(),
        ],
        critical_functional_domains=[functional_domain],
    )
    assert isinstance(f, CategoricalFusion)

    # catch and pass on validation errors
    with pytest.raises(FUSORParametersException) as excinfo:
        f = fusor_instance.fusion(
            fusion_type="CategoricalFusion", structure=[linker_element]
        )
    msg = "Fusions must contain >= 2 structural elements, or >=1 structural element and a regulatory element"
    assert msg in str(excinfo.value)

    expected = copy.deepcopy(transcript_segment_element)
    expected.exonStart = None
    with pytest.raises(FUSORParametersException) as excinfo:
        f = fusor_instance.fusion(
            structure=[
                linker_element,
                expected,
            ],
            causative_event=causative_event,
            assay=assay,
        )
    msg = "First structural element cannot be LinkerSequence"
    assert msg in str(excinfo.value)
    msg = "3' fusion partner junction must include " "starting position"
    assert msg in str(excinfo.value)

    # catch multiple errors from different validators
    with pytest.raises(FUSORParametersException) as excinfo:
        f = fusor_instance.fusion(
            structure=[
                linker_element,
                expected,
            ],
            reading_frame_preserved="not a bool",
            causative_event="other type",
        )
    msg = "Input should be a valid boolean\nInput should be a valid dictionary or instance of CausativeEvent"
    assert msg in str(excinfo.value)


@pytest.mark.asyncio()
async def test_transcript_segment_element(
    fusor_instance, transcript_segment_element, mane_transcript_segment_element
):
    """Test that transcript_segment_element method works correctly"""
    # Transcript Input
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3", exon_start=1, exon_end=8, tx_to_genomic_coords=True
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == transcript_segment_element.model_dump()

    # Genomic input, inter-residue
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        seg_start_genomic=154192135,
        seg_end_genomic=154170399,
        genomic_ac="NC_000001.11",
        tx_to_genomic_coords=False,
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == transcript_segment_element.model_dump()

    # Transcript Input
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        exon_start=1,
        exon_end=8,
        gene="TPM3",
        tx_to_genomic_coords=True,
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == transcript_segment_element.model_dump()

    expected = copy.deepcopy(transcript_segment_element)
    expected.elementGenomicStart.sequenceReference.refgetAccession = (
        "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"
    )
    expected.elementGenomicEnd.sequenceReference.refgetAccession = (
        expected.elementGenomicStart.sequenceReference.refgetAccession
    )
    expected.elementGenomicEnd.sequenceReference.id = (
        "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"
    )
    expected.elementGenomicStart.sequenceReference.id = (
        expected.elementGenomicEnd.sequenceReference.id
    )

    # Transcript Input
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        exon_start=1,
        exon_end=8,
        tx_to_genomic_coords=True,
        seq_id_target_namespace="ga4gh",
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == expected.model_dump()

    # Genomic input
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        seg_start_genomic=154192135,
        seg_end_genomic=154170399,
        genomic_ac="NC_000001.11",
        tx_to_genomic_coords=False,
        seq_id_target_namespace="ga4gh",
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == expected.model_dump()

    expected.exonEndOffset = -5
    expected.elementGenomicEnd.start = 154170404
    expected.elementGenomicEnd.end = None
    expected.elementGenomicEnd.id = "ga4gh:SL.f2Ocn2oc7X6i9fxQrRdMonLXm-W6nyn6"
    expected.elementGenomicEnd.digest = "f2Ocn2oc7X6i9fxQrRdMonLXm-W6nyn6"

    # Transcript Input
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        exon_start=1,
        exon_end=8,
        exon_end_offset=-5,
        tx_to_genomic_coords=True,
        seq_id_target_namespace="ga4gh",
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == expected.model_dump()

    # Genomic Input
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        seg_start_genomic=154192135,
        seg_end_genomic=154170404,
        genomic_ac="NC_000001.11",
        tx_to_genomic_coords=False,
        seq_id_target_namespace="ga4gh",
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == expected.model_dump()

    expected.exonEnd = None
    expected.exonEndOffset = None
    expected.elementGenomicEnd = None

    # Transcript Input
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        exon_start=1,
        tx_to_genomic_coords=True,
        seq_id_target_namespace="ga4gh",
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == expected.model_dump()

    # Genomic Input
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        seg_start_genomic=154192135,
        genomic_ac="NC_000001.11",
        tx_to_genomic_coords=False,
        seq_id_target_namespace="ga4gh",
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == expected.model_dump()

    # MANE
    tsg = await fusor_instance.transcript_segment_element(
        tx_to_genomic_coords=False,
        genomic_ac="NC_000011.10",
        seg_start_genomic=9575887,
        gene="WEE1",
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == mane_transcript_segment_element.model_dump()


def test_gene_element(fusor_instance, braf_gene_obj_min, braf_gene_obj):
    """Test that gene_element works correctly."""
    gc = fusor_instance.gene_element("BRAF", use_minimal_gene=True)
    assert gc[0]
    assert gc[1] is None
    assert isinstance(gc[0], GeneElement)
    compare_gene_obj(gc[0].gene.model_dump(), braf_gene_obj_min.model_dump())

    gc = fusor_instance.gene_element("BRAF", use_minimal_gene=False)
    assert gc[0]
    assert gc[1] is None
    assert isinstance(gc[0], GeneElement)
    compare_gene_obj(gc[0].gene.model_dump(), braf_gene_obj.model_dump())

    gc = fusor_instance.gene_element("BRA F", use_minimal_gene=True)
    assert gc[0] is None
    assert gc[1] == "gene-normalizer unable to normalize BRA F"


def test_templated_sequence_element(
    fusor_instance,
    templated_sequence_element,
):
    """Test that templated sequence element works correctly"""
    tsg = fusor_instance.templated_sequence_element(
        100, 150, "NC_000001.11", Strand.POSITIVE, coordinate_type="residue"
    )
    assert tsg.model_dump() == templated_sequence_element.model_dump()

    tsg = fusor_instance.templated_sequence_element(
        99, 150, "NC_000001.11", Strand.POSITIVE, coordinate_type="inter-residue"
    )
    assert tsg.model_dump() == templated_sequence_element.model_dump()

    # test properly defaults coordinate type to inter-residue
    tsg = fusor_instance.templated_sequence_element(
        99,
        150,
        "NC_000001.11",
        Strand.POSITIVE,
    )
    assert tsg.model_dump() == templated_sequence_element.model_dump()

    # test in-house/bespoke sequence ID
    # can't coerce namespace or translate to ga4gh ID
    with pytest.raises(IDTranslationException):
        fusor_instance.templated_sequence_element(
            200,
            300,
            "custom_ID__1",
            Strand.POSITIVE,
            coordinate_type="inter-residue",
            seq_id_target_namespace="ga4gh",
        )


def test_linker_element(fusor_instance, linker_element):
    """Test that linker_element method works correctly."""
    lc = fusor_instance.linker_element("act")
    assert lc[0]
    assert lc[1] is None
    assert lc[0].model_dump() == linker_element.model_dump()

    lc = fusor_instance.linker_element("bob!")
    assert lc[0] is None
    assert "String should match pattern '^[A-Z*\\-]*$'" in lc[1]


def test_unknown_gene_element(fusor_instance):
    """Test that unknown_gene_element method works correctly."""
    unknown_gc = fusor_instance.unknown_gene_element()
    assert unknown_gc.model_dump() == UnknownGeneElement().model_dump()


def test_multiple_possible_genes_element(fusor_instance):
    """Test that test_multiple_possible_genes_element method works correctly."""
    mult_gene = fusor_instance.multiple_possible_genes_element()
    assert mult_gene.model_dump() == MultiplePossibleGenesElement().model_dump()


def test_functional_domain(
    fusor_instance, functional_domain, functional_domain_min, functional_domain_seq_id
):
    """Test that functional_domain method works correctly"""

    def compare_domains(actual, expected):
        """Compare actual and expected functional domain data"""
        assert actual[0]
        assert actual[1] is None
        actual = actual[0].model_dump()
        expected = expected.model_dump()
        assert actual.keys() == expected.keys()
        for key in expected:
            if key == "associatedGene":
                compare_gene_obj(actual[key], expected[key])
            elif key == "sequenceLocation":
                act_sl = actual["sequenceLocation"]
                exp_sl = expected["sequenceLocation"]
                assert act_sl["id"] == exp_sl["id"]
                assert act_sl["type"] == exp_sl["type"]
                assert (
                    act_sl["sequenceReference"]["type"]
                    == exp_sl["sequenceReference"]["type"]
                )
                assert (
                    act_sl["sequenceReference"]["id"]
                    == exp_sl["sequenceReference"]["id"]
                )
                assert exp_sl.get("start") == act_sl.get("start")
                assert exp_sl.get("end") == act_sl.get("end")
            else:
                assert actual[key] == expected[key]

    cd = fusor_instance.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245",
        "BRAF",
        "NP_004324.2",
        458,
        712,
        use_minimal_gene=False,
    )
    compare_domains(cd, functional_domain)

    cd = fusor_instance.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245",
        "BRAF",
        "NP_004324.2",
        458,
        712,
        use_minimal_gene=True,
    )
    compare_domains(cd, functional_domain_min)

    cd = fusor_instance.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245",
        "BRAF",
        "NP_004324.2",
        458,
        712,
        seq_id_target_namespace="ga4gh",
        use_minimal_gene=True,
    )
    compare_domains(cd, functional_domain_seq_id)

    cd = fusor_instance.functional_domain(
        "preserveded",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245",
        "BRAF",
        "NP_004324.2",
        458,
        712,
        seq_id_target_namespace="ga4gh",
        use_minimal_gene=True,
    )
    assert cd[0] is None
    assert "Input should be 'lost' or 'preserved'" in cd[1]

    # check for protein accession
    cd = fusor_instance.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245",
        "BRAF",
        "NM_004333.4",
        458,
        712,
        seq_id_target_namespace="ga4gh",
        use_minimal_gene=True,
    )
    assert cd[0] is None
    assert "Sequence_id must be a protein accession." in cd[1]

    # check for recognized protein accession
    accession = "NP_9999.999"
    cd = fusor_instance.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245",
        "BRAF",
        accession,
        458,
        712,
        seq_id_target_namespace="ga4gh",
        use_minimal_gene=True,
    )
    assert cd[0] is None
    assert f"Accession, {accession}, not found in SeqRepo" in cd[1]

    # check that coordinates exist on sequence
    cd = fusor_instance.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245",
        "BRAF",
        "NP_004324.2",
        458,
        712000,
        seq_id_target_namespace="ga4gh",
        use_minimal_gene=True,
    )
    assert cd[0] is None
    assert (
        "End inter-residue coordinate (712000) is out of index on "
        "NP_004324.2" in cd[1]
    )


def test_regulatory_element(fusor_instance, regulatory_element, regulatory_element_min):
    """Test regulatory_element method."""

    def compare_re(actual, expected):
        """Compare actual and expected regulatory element results."""
        assert actual[0]
        assert actual[1] is None
        actual = actual[0].model_dump()
        expected = expected.model_dump()
        assert actual.keys() == expected.keys()
        assert actual["type"] == expected["type"]
        compare_gene_obj(actual["associatedGene"], expected["associatedGene"])

    re = fusor_instance.regulatory_element(RegulatoryClass.PROMOTER, "BRAF")
    compare_re(re, regulatory_element_min)

    re = fusor_instance.regulatory_element(RegulatoryClass.PROMOTER, "BRAF", False)
    compare_re(re, regulatory_element)
