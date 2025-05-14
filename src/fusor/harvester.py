"""Harvester methods for output from different fusion callers"""

import csv
from abc import ABC
from itertools import dropwhile
from pathlib import Path
from typing import ClassVar, TextIO

from civicpy.civic import FusionVariant
from pydantic import BaseModel, ConfigDict

from fusor.fusion_caller_models import (
    CIVIC,
    JAFFA,
    Arriba,
    Cicero,
    EnFusion,
    FusionCaller,
    FusionCatcher,
    Genie,
    STARFusion,
)


class FusionCallerHarvester(ABC):
    """ABC for fusion caller harvesters"""

    fusion_caller: FusionCaller
    column_rename: dict
    delimeter: str

    def _get_records(self, fusions_file: TextIO) -> csv.DictReader:
        """Read in all records from a fusions file

        :param fusion_file: The open fusions file
        :return A csv.DictReader object containing the detected fusions
        """
        return csv.DictReader(fusions_file, delimiter=self.delimeter)

    def load_records(
        self,
        fusion_path: Path,
    ) -> list[FusionCaller]:
        """Convert rows of fusion caller output to Pydantic classes

        :param fusion_path: The path to the fusions file
        :param column_rename: A dictionary of column mappings
        :param delimeter: The delimeter for the fusions file
        :raise ValueError: if the file does not exist at the specified path
        :return: A list of fusions, represented as Pydantic objects
        """
        if not fusion_path.exists():
            statement = f"{fusion_path!s} does not exist"
            raise ValueError(statement)
        fusions_list = []
        fields_to_keep = self.fusion_caller.__annotations__
        reader = self._get_records(fusion_path.open())

        for row in reader:
            filtered_row = {}
            for key, value in row.items():
                renamed_key = self.column_rename.get(key, key)
                if renamed_key in fields_to_keep:
                    filtered_row[renamed_key] = value
            fusions_list.append(self.fusion_caller(**filtered_row))
        return fusions_list


class JAFFAHarvester(FusionCallerHarvester):
    """Class for harvesting JAFFA data"""

    column_rename: ClassVar[dict[str, str]] = {
        "fusion genes": "fusion_genes",
        "spanning reads": "spanning_reads",
        "spanning pairs": "spanning_pairs",
    }
    delimeter = ","
    fusion_caller = JAFFA


class StarFusionHarvester(FusionCallerHarvester):
    """Class for harvesting STAR-Fusion data"""

    column_rename: ClassVar[dict[str, str]] = {
        "LeftGene": "left_gene",
        "RightGene": "right_gene",
        "LeftBreakpoint": "left_breakpoint",
        "RightBreakpoint": "right_breakpoint",
        "JunctionReadCount": "junction_read_count",
        "SpanningFragCount": "spanning_frag_count",
    }
    delimeter = "\t"
    fusion_caller = STARFusion


class FusionCatcherHarvester(FusionCallerHarvester):
    """Class for harvesting FusionCatcher data"""

    column_rename: ClassVar[dict[str, str]] = {
        "Gene_1_symbol(5end_fusion_partner)": "five_prime_partner",
        "Gene_2_symbol(3end_fusion_partner)": "three_prime_partner",
        "Fusion_point_for_gene_1(5end_fusion_partner)": "five_prime_fusion_point",
        "Fusion_point_for_gene_2(3end_fusion_partner)": "three_prime_fusion_point",
        "Predicted_effect": "predicted_effect",
        "Spanning_unique_reads": "spanning_unique_reads",
        "Spanning_pairs": "spanning_reads",
        "Fusion_sequence": "fusion_sequence",
    }
    delimeter = "\t"
    fusion_caller = FusionCatcher


class ArribaHarvester(FusionCallerHarvester):
    """Class for harvesting Arriba data"""

    column_rename: ClassVar[dict[str, str]] = {
        "#gene1": "gene1",
        "strand1(gene/fusion)": "strand1",
        "strand2(gene/fusion)": "strand2",
        "type": "event_type",
        "reading_frame": "rf",
    }
    delimeter = "\t"
    fusion_caller = Arriba


class CiceroHarvester(FusionCallerHarvester):
    """Class for harvesting Cicero data"""

    column_rename: ClassVar[dict[str, str]] = {
        "geneA": "gene_5prime",
        "geneB": "gene_3prime",
        "chrA": "chr_5prime",
        "chrB": "chr_3prime",
        "posA": "pos_5prime",
        "posB": "pos_3prime",
        "type": "event_type",
        "readsA": "reads_5prime",
        "readsB": "reads_3prime",
        "coverageA": "coverage_5prime",
        "coverageB": "coverage_3prime",
    }
    delimeter = "\t"
    fusion_caller = Cicero


class EnFusionHarvester(FusionCallerHarvester):
    """Class for harvesting EnFusion data"""

    column_rename: ClassVar[dict[str, str]] = {
        "Gene1": "gene_5prime",
        "Gene2": "gene_3prime",
        "Chr1": "chr_5prime",
        "Chr2": "chr_3prime",
        "Break1": "break_5prime",
        "Break2": "break_3prime",
        "FusionJunctionSequence": "fusion_junction_sequence",
    }
    delimeter = "\t"
    fusion_caller = EnFusion

    def _get_records(self, fusions_file: TextIO) -> csv.DictReader:
        """Read in all records from a fusions file

        :param fusion_file: The open fusions file
        :return A csv.DictReader object containing the detected fusions
        """
        fusion_lines = dropwhile(
            lambda line: not line.startswith("UnorderedFusion"), fusions_file
        )
        return csv.DictReader(fusion_lines, delimiter=self.delimeter)


class GenieHarvester(FusionCallerHarvester):
    """Class for harvesting Genie data"""

    column_rename: ClassVar[dict[str, str]] = {
        "Site1_Hugo_Symbol": "site1_hugo",
        "Site2_Hugo_Symbol": "site2_hugo",
        "Site1_Chromosome": "site1_chrom",
        "Site2_Chromosome": "site2_chrom",
        "Site1_Position": "site1_pos",
        "Site2_Position": "site2_pos",
        "Site2_Effect_On_Frame": "reading_frame",
        "Annotation": "annot",
    }
    delimeter = "\t"
    fusion_caller = Genie


class CIVICHarvester(BaseModel):
    """Class for harvesting CIViC Fusion objects"""

    fusions_list: list[FusionVariant]
    model_config = ConfigDict(arbitrary_types_allowed=True)

    def load_records(self) -> list[CIVIC]:
        """Extract data from CIVIC fusion objects

        :return A list of CIVIC objects
        """
        processed_fusions = []
        for fusion in self.fusions_list:
            params = {
                "vicc_compliant_name": fusion.vicc_compliant_name,
                "five_prime_end_exon_coords": fusion.five_prime_end_exon_coordinates,
                "three_prime_start_exon_coords": fusion.three_prime_start_exon_coordinates,
                "molecular_profiles": fusion.molecular_profiles,
            }
            processed_fusions.append(CIVIC(**params))
        return processed_fusions
