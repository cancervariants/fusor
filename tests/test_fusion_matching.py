"""Module for testing fusion matching module"""

from pathlib import Path

import pytest
from civicpy import civic
from cool_seq_tool.schemas import Assembly, CoordinateType

from fusor.fusion_matching import FusionSources
from fusor.harvester import CIVICHarvester, StarFusionHarvester


@pytest.mark.asyncio()
async def test_fusion_matching_valid(
    fixture_data_dir, translator_instance, fusion_matching_instance
):
    """Test fusion matching worklow using example output from STAR-Fusion"""
    path = Path(fixture_data_dir / "star-fusion.fusion_predictions.abridged.tsv")
    harvester = StarFusionHarvester()
    fusions_list = harvester.load_records(path)

    # Test KIF5B::RET Fusion
    assayed_fusion_star_fusion = await translator_instance.from_star_fusion(
        fusions_list[0], CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    civic_variants = civic.get_all_fusion_variants(include_status="accepted")
    harvester = CIVICHarvester()
    harvester.fusions_list = civic_variants
    fusions = harvester.load_records()

    matching_categorical_fusions = [
        await translator_instance.from_civic(fusion)
        for fusion in fusions
        if fusion.vicc_compliant_name
        in {
            "KIF5B(entrez:3799)::RET(entrez:5979)",
            "v::RET(entrez:5979)",
        }
    ]
    matches = await fusion_matching_instance.match_fusion(
        assayed_fusion=assayed_fusion_star_fusion,
        sources=[FusionSources.CIVIC],
        cache_dir=Path("src/fusor/data/"),
    )
    assert len(matches) == 2
    assert matches[0][0] == matching_categorical_fusions[0]
    assert matches[0][1] == 10
    assert matches[1][0] == matching_categorical_fusions[1]
    assert matches[1][1] == 1

    fusion_matching_instance.categorical_fusions.clear()

    # Test EML4::ALK fusion (exon 13 of EML4 fused with exon 20 of ALK)
    # Note: There are numerous EML4::ALK fusions in CIViC, but only one entry that
    # contains the exact matching transcript junctions
    assayed_fusion_star_fusion = await translator_instance.from_star_fusion(
        fusions_list[1], CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    matching_categorical_fusions = [
        await translator_instance.from_civic(fusion)
        for fusion in fusions
        if fusion.vicc_compliant_name
        in {
            "EML4(entrez:27436)::ALK(entrez:238)",
            "v::ALK(entrez:238)",
        }
    ]
    matches = await fusion_matching_instance.match_fusion(
        assayed_fusion=assayed_fusion_star_fusion,
        sources=[FusionSources.CIVIC],
        cache_dir=Path("src/fusor/data/"),
    )
    assert len(matches) == 2
    assert matches[0][0] == matching_categorical_fusions[0]
    assert matches[0][1] == 10
    assert matches[1][0] == matching_categorical_fusions[1]
    assert matches[1][1] == 5


@pytest.mark.asyncio()
async def test_fusion_matching_invalid(
    fixture_data_dir, translator_instance, fusion_matching_instance
):
    """Test case where assayed fusion would not have a categorical fusion match"""
    path = Path(fixture_data_dir / "star-fusion.fusion_predictions.abridged.tsv")
    harvester = StarFusionHarvester()
    fusions_list = harvester.load_records(path)

    assayed_fusion_star_fusion = await translator_instance.from_star_fusion(
        fusions_list[36], CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    matches = await fusion_matching_instance.match_fusion(
        assayed_fusion=assayed_fusion_star_fusion,
        sources=[FusionSources.CIVIC],
        cache_dir=Path("src/fusor/data/"),
    )
    assert not matches
