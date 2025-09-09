"""Module for testing fusion caller harvesters"""

from pathlib import Path

import pytest
from cool_seq_tool.schemas import Assembly

from fusor.harvester import (
    ArribaHarvester,
    CiceroHarvester,
    CIVICHarvester,
    EnFusionHarvester,
    FusionCatcherHarvester,
    GenieHarvester,
    JAFFAHarvester,
    MOAHarvester,
    StarFusionHarvester,
)


async def test_get_jaffa_records(fixture_data_dir, fusor_instance):
    """Test that JAFFAHarvester works correctly"""
    path = Path(fixture_data_dir / "jaffa_results_test.csv")
    harvester = JAFFAHarvester(fusor_instance, assembly=Assembly.GRCH38)
    records = await harvester.load_records(path)
    assert len(records) == 2

    records = await harvester.load_record_table(path)
    assert len(records) == 2

    path = Path(fixture_data_dir / "jaffa_resultss_test.csv")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert await harvester.load_records(path)


async def test_get_star_fusion_records(fixture_data_dir, fusor_instance):
    """Test that STARFusionHarvester works correctly"""
    path = Path(fixture_data_dir / "star_fusion_test.tsv")
    harvester = StarFusionHarvester(fusor_instance, assembly=Assembly.GRCH38)
    records = await harvester.load_records(path)
    assert len(records) == 3

    records = await harvester.load_record_table(path)
    assert len(records) == 3

    path = Path(fixture_data_dir / "star_fusion_test.tsvs")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert await harvester.load_records(path)


async def test_get_fusion_catcher_records(fixture_data_dir, fusor_instance):
    """Test that FusionCatcherHarvester works correctly"""
    path = Path(fixture_data_dir / "fusion_catcher_test.txt")
    harvester = FusionCatcherHarvester(fusor_instance, assembly=Assembly.GRCH38)
    fusions_list = await harvester.load_records(path)
    assert len(fusions_list) == 3

    records = await harvester.load_record_table(path)
    assert len(records) == 3

    path = Path(fixture_data_dir / "fusionn_catcher.txts")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert await harvester.load_records(path)


async def test_get_arriba_records(fixture_data_dir, fusor_instance):
    """Test that ArribaHarvester works correctly"""
    path = Path(fixture_data_dir / "fusions_arriba_test.tsv")
    harvester = ArribaHarvester(fusor_instance, assembly=Assembly.GRCH37)
    fusions_list = await harvester.load_records(path)
    assert len(fusions_list) == 1

    records = await harvester.load_record_table(path)
    assert len(records) == 1

    # Test ITD example
    path = Path(fixture_data_dir / "arriba_fusions_itd.tsv")
    harvester = ArribaHarvester(fusor_instance, assembly=Assembly.GRCH38)
    itd_list = await harvester.load_record_table(path)
    assert len(itd_list) == 1

    path = Path(fixture_data_dir / "fusionss_arriba_test.tsv")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert await harvester.load_records(path)


async def test_get_cicero_records(fixture_data_dir, fusor_instance):
    """Test that CiceroHarvester works correctly"""
    path = Path(fixture_data_dir / "annotated.fusion.txt")
    harvester = CiceroHarvester(fusor_instance, assembly=Assembly.GRCH38)
    fusions_list = await harvester.load_records(path)
    assert len(fusions_list) == 1

    records = await harvester.load_record_table(path)
    assert len(records) == 1

    path = Path(fixture_data_dir / "annnotated.fusion.txt")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert await harvester.load_records(path)


async def test_get_enfusion_records(fixture_data_dir, fusor_instance):
    """Test that EnFusionHarvester works correctly"""
    path = Path(fixture_data_dir / "enfusion_test.csv")
    harvester = EnFusionHarvester(fusor_instance, assembly=Assembly.GRCH38)
    fusions_list = await harvester.load_records(path)
    assert len(fusions_list) == 1

    records = await harvester.load_record_table(path)
    assert len(records) == 1

    path = Path(fixture_data_dir / "enfusions_test.csv")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert await harvester.load_records(path)


async def test_get_genie_records(fixture_data_dir, fusor_instance):
    """Test that GenieHarvester works correctly"""
    path = Path(fixture_data_dir / "genie_test.txt")
    harvester = GenieHarvester(fusor_instance, assembly=Assembly.GRCH38)
    fusions_list = await harvester.load_records(path)
    assert len(fusions_list) == 1

    records = await harvester.load_record_table(path)
    assert len(records) == 1

    path = Path(fixture_data_dir / "genie_tests.txt")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert await harvester.load_records(path)


async def test_get_civic_records(fixture_data_dir, fusor_instance):
    """Test that CIVICHarvester works correctly"""
    harvester = CIVICHarvester(
        fusor_instance, local_cache_path=fixture_data_dir / "civic_cache.pkl"
    )
    harvester.fusions_list = harvester.fusions_list[
        :5
    ]  # Look at first 5 records in test
    fusions_list = await harvester.load_records()
    assert len(fusions_list) == 5


def test_get_moa_records(fusor_instance, fixture_data_dir):
    """Test that MOAHarvester works correctly"""
    harvester = MOAHarvester(fusor_instance, cache_dir=fixture_data_dir, use_local=True)
    fusions_list = harvester.load_records()
    assert len(fusions_list) == 97
