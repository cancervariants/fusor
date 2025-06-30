"""Module for testing fusion caller classes"""

from pathlib import Path

import pytest
from civicpy import civic
from cool_seq_tool.schemas import Assembly

from fusor.harvester import (
    ArribaHarvester,
    CiceroHarvester,
    CIVICHarvester,
    EnFusionHarvester,
    FusionCatcherHarvester,
    GenieHarvester,
    JAFFAHarvester,
    StarFusionHarvester,
)


async def test_get_jaffa_records(fixture_data_dir, translator_instance):
    """Test that get_jaffa_records works correctly"""
    path = Path(fixture_data_dir / "jaffa_results_test.csv")
    harvester = JAFFAHarvester(translator_instance, assembly=Assembly.GRCH38.value)
    records = await harvester.load_records(path)
    assert len(records) == 9

    path = Path(fixture_data_dir / "jaffa_resultss_test.csv")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert await harvester.load_records(path)


async def test_get_star_fusion_records(fixture_data_dir, translator_instance):
    """Test that get_star_fusion_records works correctly"""
    path = Path(fixture_data_dir / "star-fusion.fusion_predictions.abridged.tsv")
    harvester = StarFusionHarvester(translator_instance, assembly=Assembly.GRCH38.value)
    records = await harvester.load_records(path)
    assert len(records) == 37

    path = Path(fixture_data_dir / "star-fusion.fusion_predictions.abridged.tsvs")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert await harvester.load_records(path)


async def test_get_fusion_catcher_records(fixture_data_dir, translator_instance):
    """Test that get_fusion_catcher_records works correctly"""
    path = Path(fixture_data_dir / "fusion_catcher_test.txt")
    harvester = FusionCatcherHarvester(
        translator_instance, assembly=Assembly.GRCH38.value
    )
    fusions_list = await harvester.load_records(path)
    assert len(fusions_list) == 14

    path = Path(fixture_data_dir / "fusionn_catcher.txts")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert await harvester.load_records(path)


async def test_get_arriba_records(fixture_data_dir, translator_instance):
    """Test that get_arriba_records works correctly"""
    path = Path(fixture_data_dir / "fusions_arriba_test.tsv")
    harvester = ArribaHarvester(translator_instance, assembly=Assembly.GRCH37.value)
    fusions_list = await harvester.load_records(path)
    assert len(fusions_list) == 1

    path = Path(fixture_data_dir / "fusionss_arriba_test.tsv")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert await harvester.load_records(path)


async def test_get_cicero_records(fixture_data_dir, translator_instance):
    """Test that get_cicero_records works correctly"""
    path = Path(fixture_data_dir / "annotated.fusion.txt")
    harvester = CiceroHarvester(translator_instance, assembly=Assembly.GRCH38.value)
    fusions_list = await harvester.load_records(path)
    assert len(fusions_list) == 1

    path = Path(fixture_data_dir / "annnotated.fusion.txt")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert await harvester.load_records(path)


async def test_get_enfusion_records(fixture_data_dir, translator_instance):
    """Test that get_enfusion_records works correctly"""
    path = Path(fixture_data_dir / "enfusion_test.csv")
    harvester = EnFusionHarvester(translator_instance, assembly=Assembly.GRCH38.value)
    fusions_list = await harvester.load_records(path)
    assert len(fusions_list) == 1

    path = Path(fixture_data_dir / "enfusions_test.csv")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert await harvester.load_records(path)


async def test_get_genie_records(fixture_data_dir, translator_instance):
    """Test that get_genie_records works correctly"""
    path = Path(fixture_data_dir / "genie_test.txt")
    harvester = GenieHarvester(translator_instance, assembly=Assembly.GRCH38.value)
    fusions_list = await harvester.load_records(path)
    assert len(fusions_list) == 1

    path = Path(fixture_data_dir / "genie_tests.txt")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert await harvester.load_records(path)


async def test_get_civic_records(translator_instance):
    """Test that get_civic_records works correctly"""
    civic_variants = civic.get_all_fusion_variants()
    civic_variants = civic_variants[:5]  # Look at first 5 records in test
    harvester = CIVICHarvester(translator=translator_instance)
    harvester.fusions_list = civic_variants
    fusions_list = await harvester.load_records()
    assert len(fusions_list) == 5
