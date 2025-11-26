"""Module for testing build_fusion method"""

import pytest
from cool_seq_tool.schemas import AnnotationLayer, ManeStatus

from fusor.models import AssayedFusion, CategoricalFusion
from fusor.transcript_junction_builder import TranscriptJunctionBuilder


def fusion_example_bcr_abl1(**kwargs):
    """Create example `AssayedFusion` object for BCR::ABL1"""
    params = {
        "type": "AssayedFusion",
        "structure": [
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_004327.3",
                "transcriptStatus": "longest_compatible_remaining",
                "strand": 1,
                "exonEnd": 16,
                "exonEndOffset": -12,
                "gene": {
                    "primaryCoding": {
                        "id": "hgnc:1014",
                        "code": "HGNC:1014",
                        "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                    },
                    "conceptType": "Gene",
                    "name": "BCR",
                },
                "elementGenomicEnd": {
                    "id": "ga4gh:SL.tQkBo8Xr_nDIElZFWqsy44yXXLxjW5Rd",
                    "type": "SequenceLocation",
                    "digest": "tQkBo8Xr_nDIElZFWqsy44yXXLxjW5Rd",
                    "sequenceReference": {
                        "id": "refseq:NC_000022.11",
                        "type": "SequenceReference",
                        "refgetAccession": "SQ.7B7SHsmchAR0dFcDCuSFjJAo7tX87krQ",
                    },
                    "end": 23295143,
                    "extensions": [{"name": "is_exonic", "value": True}],
                },
            },
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_005157.5",
                "transcriptStatus": "longest_compatible_remaining",
                "strand": 1,
                "exonStart": 11,
                "exonStartOffset": 321,
                "gene": {
                    "primaryCoding": {
                        "id": "hgnc:76",
                        "code": "HGNC:76",
                        "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                    },
                    "conceptType": "Gene",
                    "name": "ABL1",
                },
                "elementGenomicStart": {
                    "id": "ga4gh:SL.4AalyLSqDQxfxmlLmx0T9dizutqExADN",
                    "type": "SequenceLocation",
                    "digest": "4AalyLSqDQxfxmlLmx0T9dizutqExADN",
                    "sequenceReference": {
                        "id": "refseq:NC_000009.12",
                        "type": "SequenceReference",
                        "refgetAccession": "SQ.KEO-4XBcm1cxeo_DIQ8_ofqGUkp4iZhI",
                    },
                    "start": 130884289,
                    "extensions": [{"name": "is_exonic", "value": True}],
                },
            },
        ],
    }
    assayed_fusion = AssayedFusion(**params)
    return assayed_fusion.model_copy(update=kwargs)


def fusion_example_bcr_abl1_gene_elements(**kwargs):
    """Create example `AssayedFusion` object for BCR::ABL1 using only GeneElements"""
    params = {
        "type": "AssayedFusion",
        "structure": [
            {
                "type": "GeneElement",
                "gene": {
                    "primaryCoding": {
                        "id": "hgnc:1014",
                        "code": "HGNC:1014",
                        "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                    },
                    "conceptType": "Gene",
                    "name": "BCR",
                },
            },
            {
                "type": "GeneElement",
                "gene": {
                    "primaryCoding": {
                        "id": "hgnc:76",
                        "code": "HGNC:76",
                        "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                    },
                    "conceptType": "Gene",
                    "name": "ABL1",
                },
            },
        ],
    }
    assayed_fusion = AssayedFusion(**params)
    return assayed_fusion.model_copy(update=kwargs)


def fusion_example_eml4_alk(**kwargs):
    """Create example `CategoricalFusion` object for EML4::ALK"""
    params = {
        "type": "CategoricalFusion",
        "structure": [
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_019063.4",
                "transcriptStatus": "longest_compatible_remaining",
                "strand": 1,
                "exonEnd": 20,
                "exonEndOffset": 0,
                "gene": {
                    "primaryCoding": {
                        "id": "hgnc:1316",
                        "code": "HGNC:1316",
                        "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                    },
                    "conceptType": "Gene",
                    "name": "EML4",
                },
                "elementGenomicEnd": {
                    "id": "ga4gh:SL.-IgV899bqi2cN3ugOPhJG_ZAbuUgrN7N",
                    "type": "SequenceLocation",
                    "digest": "-IgV899bqi2cN3ugOPhJG_ZAbuUgrN7N",
                    "sequenceReference": {
                        "id": "refseq:NC_000002.12",
                        "type": "SequenceReference",
                        "refgetAccession": "SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
                    },
                    "end": 42325554,
                    "extensions": [{"name": "is_exonic", "value": True}],
                },
            },
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_004304.4",
                "transcriptStatus": "longest_compatible_remaining",
                "strand": -1,
                "exonStart": 20,
                "exonStartOffset": 0,
                "gene": {
                    "primaryCoding": {
                        "id": "hgnc:427",
                        "code": "HGNC:427",
                        "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                    },
                    "conceptType": "Gene",
                    "name": "ALK",
                },
                "elementGenomicStart": {
                    "id": "ga4gh:SL.Eu_igVd9zOahn3tFN-pyxtphUmrSlRAh",
                    "type": "SequenceLocation",
                    "digest": "Eu_igVd9zOahn3tFN-pyxtphUmrSlRAh",
                    "sequenceReference": {
                        "id": "refseq:NC_000002.12",
                        "type": "SequenceReference",
                        "refgetAccession": "SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
                    },
                    "end": 29223528,
                    "extensions": [{"name": "is_exonic", "value": True}],
                },
            },
        ],
    }
    categorical_fusion = CategoricalFusion(**params)
    return categorical_fusion.model_copy(update=kwargs)


test_data = [
    (
        "BCR",
        "ABL1",
        AnnotationLayer.CDNA,
        "NM_004327.3",
        "NM_005157.5",
        3000,
        2000,
        True,
        fusion_example_bcr_abl1(),
    ),
    (
        "BCR",
        "ABL1",
        AnnotationLayer.GENOMIC,
        "NC_000022.11",
        "NC_000009.12",
        23295143,
        130884290,
        True,
        fusion_example_bcr_abl1(),
    ),
    (
        "BCR",
        "ABL1",
        None,
        None,
        None,
        None,
        None,
        True,
        fusion_example_bcr_abl1_gene_elements(),
    ),
    (
        "EML4",
        "ALK",
        AnnotationLayer.CDNA,
        "NM_019063.4",
        "NM_004304.4",
        2242,
        3173,
        False,
        fusion_example_eml4_alk(),
    ),
    (
        "EML4",
        "ALK",
        AnnotationLayer.GENOMIC,
        "NC_000002.12",
        "NC_000002.12",
        42325554,
        29223528,
        False,
        fusion_example_eml4_alk(),
    ),
]


def increment_refseq(refseq: str) -> str:
    """Increment RefSeq transcript version for genomic tests

    :param refseq: A RefSeq accession
    :return: An incremented RefSeq accession
    """
    base, version = refseq.split(".")
    version = int(version) + 1
    return base + "." + str(version)


@pytest.mark.asyncio
@pytest.mark.parametrize(
    (
        "five_prime_gene",
        "three_prime_gene",
        "annotation_layer",
        "five_prime_reference_sequence",
        "three_prime_reference_sequence",
        "five_prime_junction",
        "three_prime_junction",
        "assayed_fusion",
        "expected",
    ),
    test_data,
)
async def test_build_fusion(
    fusor_instance,
    five_prime_gene,
    three_prime_gene,
    annotation_layer,
    five_prime_reference_sequence,
    three_prime_reference_sequence,
    five_prime_junction,
    three_prime_junction,
    assayed_fusion,
    expected,
):
    tbf = TranscriptJunctionBuilder(
        fusor=fusor_instance,
        five_prime_gene=five_prime_gene,
        three_prime_gene=three_prime_gene,
        five_prime_reference_sequence=five_prime_reference_sequence,
        three_prime_reference_sequence=three_prime_reference_sequence,
        annotation_type=annotation_layer,
        five_prime_junction=five_prime_junction,
        three_prime_junction=three_prime_junction,
        assayed_fusion=assayed_fusion,
    )

    fusion = await tbf.build_fusion()
    if annotation_layer == AnnotationLayer.GENOMIC:
        expected.structure[0].transcript = increment_refseq(
            expected.structure[0].transcript
        )
        expected.structure[0].transcriptStatus = ManeStatus.SELECT
        expected.structure[1].transcript = increment_refseq(
            expected.structure[1].transcript
        )
        expected.structure[1].transcriptStatus = ManeStatus.SELECT
    assert fusion.type == expected.type
    assert fusion.structure == expected.structure


test_data_invalid = [
    (
        "BCR",
        "ABL1",
        AnnotationLayer.PROTEIN,
        "NM_004327.3",
        "NM_005157.5",
        3000,
        2000,
        True,
        "Only c. or g. RefSeq accessions are supported",
        ValueError,
    ),
    (
        "BCR",
        "ABL1",
        AnnotationLayer.CDNA,
        "NM_004327.3",
        "NM_005157.5",
        None,
        2000,
        True,
        "Please provide 5' junction location",
        ValueError,
    ),
    (
        "BCR",
        "ABL1",
        AnnotationLayer.CDNA,
        "NM_004327.3",
        "NM_005157.5",
        3000,
        None,
        True,
        "Please provide 3' junction location",
        ValueError,
    ),
]


@pytest.mark.parametrize(
    (
        "five_prime_gene",
        "three_prime_gene",
        "annotation_layer",
        "five_prime_reference_sequence",
        "three_prime_reference_sequence",
        "five_prime_junction",
        "three_prime_junction",
        "assayed_fusion",
        "error_statement",
        "expected",
    ),
    test_data_invalid,
)
def test_build_fusion_invalid(
    fusor_instance,
    five_prime_gene,
    three_prime_gene,
    annotation_layer,
    five_prime_reference_sequence,
    three_prime_reference_sequence,
    five_prime_junction,
    three_prime_junction,
    assayed_fusion,
    error_statement,
    expected,
):
    """Test invalid cases"""
    with pytest.raises(
        expected,
        match=error_statement,
    ):
        TranscriptJunctionBuilder(
            fusor=fusor_instance,
            five_prime_gene=five_prime_gene,
            three_prime_gene=three_prime_gene,
            five_prime_reference_sequence=five_prime_reference_sequence,
            three_prime_reference_sequence=three_prime_reference_sequence,
            annotation_type=annotation_layer,
            five_prime_junction=five_prime_junction,
            three_prime_junction=three_prime_junction,
            assayed_fusion=assayed_fusion,
        )
