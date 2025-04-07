"""Module for testing fusion caller classes"""

from pathlib import Path

import pytest
from civicpy import civic

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


def test_get_jaffa_records(fixture_data_dir):
    """Test that get_jaffa_records works correctly"""
    path = Path(fixture_data_dir / "jaffa_results.csv")
    harvester = JAFFAHarvester()
    records = harvester.load_records(path)
    assert len(records) == 491

    path = Path(fixture_data_dir / "jaffa_resultss.csv")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert harvester.load_records(path)


def test_get_star_fusion_records(fixture_data_dir):
    """Test that get_star_fusion_records works correctly"""
    path = Path(fixture_data_dir / "star-fusion.fusion_predictions.abridged.tsv")
    harvester = StarFusionHarvester()
    records = harvester.load_records(path)
    assert len(records) == 37

    path = Path(fixture_data_dir / "star-fusion.fusion_predictions.abridged.tsvs")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert harvester.load_records(path)


def test_get_fusion_catcher_records(fixture_data_dir):
    """Test that get_fusion_catcher_records works correctly"""
    path = Path(fixture_data_dir / "final-list_candidate-fusion-genes.txt")
    harvester = FusionCatcherHarvester()
    fusions_list = harvester.load_records(path)
    assert len(fusions_list) == 355

    path = Path(fixture_data_dir / "final-list_candidate-fusion-genes.txts")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert harvester.load_records(path)


def test_get_arriba_records(fixture_data_dir):
    """Test that get_arriba_records works correctly"""
    path = Path(fixture_data_dir / "fusions_arriba_test.tsv")
    harvester = ArribaHarvester()
    fusions_list = harvester.load_records(path)
    assert len(fusions_list) == 1

    path = Path(fixture_data_dir / "fusionsd_arriba_test.tsv")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert harvester.load_records(path)


def test_get_cicero_records(fixture_data_dir):
    """Test that get_cicero_records works correctly"""
    path = Path(fixture_data_dir / "annotated.fusion.txt")
    harvester = CiceroHarvester()
    fusions_list = harvester.load_records(path)
    assert len(fusions_list) == 1

    path = Path(fixture_data_dir / "annnotated.fusion.txt")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert harvester.load_records(path)


def test_get_enfusion_records(fixture_data_dir):
    """Test that get_enfusion_records works correctly"""
    path = Path(fixture_data_dir / "enfusion_test.csv")
    harvester = EnFusionHarvester()
    fusions_list = harvester.load_records(path)
    assert len(fusions_list) == 1

    path = Path(fixture_data_dir / "enfusions_test.csv")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert harvester.load_records(path)


def test_get_genie_records(fixture_data_dir):
    """Test that get_genie_records works correctly"""
    path = Path(fixture_data_dir / "genie_test.txt")
    harvester = GenieHarvester()
    fusions_list = harvester.load_records(path)
    assert len(fusions_list) == 1

    path = Path(fixture_data_dir / "genie_tests.txt")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert harvester.load_records(path)


def test_get_civic_records():
    """Test that get_civic_records works correctly"""
    civic_variants = civic.get_all_fusion_variants()
    harvester = CIVICHarvester(fusions_list=civic_variants)
    fusions_list = harvester.load_records()
    assert len(fusions_list) == len(civic_variants)
