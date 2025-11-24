"""Module for testing quick_start method"""

import pytest

from fusor.models import AssayedFusion, CategoricalFusion
from fusor.quick_start import FUSORQuickStart


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
        "NM_004327.3:c.3000",
        "NM_005157.5:c.2000",
        "BCR",
        "ABL1",
        fusion_example_bcr_abl1(),
    ),
    (
        "NM_019063.4:c.2242",
        "NM_004304.4:c.3173",
        "EML4",
        "ALK",
        fusion_example_eml4_alk(),
    ),
]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    (
        "five_prime_junction",
        "three_prime_junction",
        "five_prime_gene",
        "three_prime_gene",
        "expected",
    ),
    test_data,
)
async def test_quick_start(
    fusor_instance,
    five_prime_junction,
    three_prime_junction,
    five_prime_gene,
    three_prime_gene,
    expected,
):
    fqs = FUSORQuickStart(
        fusor=fusor_instance,
        five_prime_junction=five_prime_junction,
        three_prime_junction=three_prime_junction,
        five_prime_gene=five_prime_gene,
        three_prime_gene=three_prime_gene,
    )

    if isinstance(expected, ValueError):
        with pytest.raises(
            expected,
            match="The fusion junction locations must be described using c. coordinates",
        ):
            await fqs.quick_start()

    fusion = await fqs.quick_start()
    assert fusion.structure == expected.structure


test_data_invalid = [
    ("NM_019063.4:c.2242", "NC_000002.12:g.29223528", "EML4", "ALK", ValueError),
]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    (
        "five_prime_junction",
        "three_prime_junction",
        "five_prime_gene",
        "three_prime_gene",
        "expected",
    ),
    test_data_invalid,
)
def test_quick_start_invali(
    fusor_instance,
    five_prime_junction,
    three_prime_junction,
    five_prime_gene,
    three_prime_gene,
    expected,
):
    with pytest.raises(
        expected,
        match="The fusion junction locations must be described using c. coordinates",
    ):
        FUSORQuickStart(
            fusor=fusor_instance,
            five_prime_junction=five_prime_junction,
            three_prime_junction=three_prime_junction,
            five_prime_gene=five_prime_gene,
            three_prime_gene=three_prime_gene,
        )
