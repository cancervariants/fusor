"""Module for testing fusion matching module"""

from pathlib import Path

import pytest
from civicpy import civic
from cool_seq_tool.schemas import Assembly, CoordinateType

from fusor.fusion_matching import CIVICCategoricalFusions
from fusor.harvester import CIVICHarvester, StarFusionHarvester


@pytest.mark.asyncio()
async def test_fusion_matching(
    fixture_data_dir, translator_instance, fusion_matching_instance
):
    """Test fusion matching worklow using example output from STAR-Fusion"""
    path = Path(fixture_data_dir / "star-fusion.fusion_predictions.abridged.tsv")
    harvester = StarFusionHarvester()
    fusions_list = harvester.load_records(path)
    assayed_fusion_star_fusion = await translator_instance.from_star_fusion(
        fusions_list[0], CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    harvester = CIVICHarvester(
        fusions_list=civic.get_all_fusion_variants(include_status="accepted")
    )
    fusions = harvester.load_records()

    # Test valid fusion that exists in CIViC
    matching_categorical_fusions = [
        await translator_instance.from_civic(fusion)
        for fusion in fusions
        if fusion.vicc_compliant_name
        in {
            "KIF5B(entrez:3799)::RET(entrez:5979)",
            "v::RET(entrez:5979)",
        }
    ]
    path = Path("src/fusor/data/civic_translated_fusions.pkl")
    matches = await fusion_matching_instance.match_fusion(
        assayed_fusion_star_fusion,
        [(CIVICCategoricalFusions(translator_instance), path)],
    )
    assert len(matches) == 2
    assert matches[0][0] == matching_categorical_fusions[0]
    assert matches[0][1] == 10
    assert matches[1][0] == matching_categorical_fusions[1]
    assert matches[1][1] == 1

    # Test case where fusion does not exist in CIViC
    assayed_fusion_star_fusion = await translator_instance.from_star_fusion(
        fusions_list[36], CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    matches = await fusion_matching_instance.match_fusion(
        assayed_fusion_star_fusion,
        [(CIVICCategoricalFusions(translator_instance), path)],
    )
    assert not matches
