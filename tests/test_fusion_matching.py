"""Module for testing fusion matching module"""

from pathlib import Path

import pytest
import yaml
from cool_seq_tool.schemas import Assembly

from fusor.fusion_matching import MatchType
from fusor.harvester import StarFusionHarvester


def _assert_subset(actual: dict, expected: dict) -> None:
    """Compare actual and expected categorical fusion fields

    :param acutal: The actual categorical fusion dictionary
    :param expected: The expected categorical fusion dictionary
    """
    for i, expected_item in enumerate(expected):
        actual_item = actual[i]
        for field in expected_item:
            expected_value = expected_item[field]
            actual_value = actual_item[field]

            if isinstance(expected_value, dict):
                for subfield in expected_value:
                    assert actual_value[subfield] == expected_value[subfield]
            else:
                assert actual_value == expected_value


@pytest.mark.asyncio
async def test_fusion_matching(
    fixture_data_dir, fusion_matching_instance, fusor_instance
):
    """Test fusion matching worklow using example output from STAR-Fusion"""
    with Path.open(
        Path(__file__).parent / "fusion_matching_test_cases.yaml"
    ) as test_cases:
        test_cases = yaml.safe_load(test_cases)["tests"]

    # Load STAR-Fusion records
    path = Path(fixture_data_dir / "star_fusion_test.tsv")
    harvester = StarFusionHarvester(fusor=fusor_instance, assembly=Assembly.GRCH38)
    fusions_list = await harvester.load_records(path)
    test_list = fusions_list[-2:]

    for case in test_cases:
        assayed_fusion = test_list[case["input_index"]]
        fusion_matching_instance.assayed_fusions = [assayed_fusion]
        matches = await fusion_matching_instance.match_fusion()

        if not case["expected_matches"]:
            assert not matches[0]
            continue

        matches = matches[0]
        for i, expected in enumerate(case["expected_matches"]):
            fusion, match_type = matches[i]
            expected_fusion = expected["fields"]
            expected_match_type = MatchType[expected["match_type"]]

            fusion = fusion.model_dump(exclude_none=True)
            for field, expected_value in expected_fusion.items():
                if field == "structure":
                    _assert_subset(fusion["structure"], expected_value)
                else:
                    assert fusion[field] == expected_value

            assert match_type == expected_match_type
