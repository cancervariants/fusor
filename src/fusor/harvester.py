"""Harvester methods for output from different fusion callers"""

import csv
import logging
from abc import ABC
from itertools import dropwhile
from pathlib import Path
from typing import ClassVar, TextIO

from civicpy import civic
from cool_seq_tool.schemas import Assembly, CoordinateType

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
from fusor.fusor import FUSOR
from fusor.models import AssayedFusion, CategoricalFusion
from fusor.translator import (
    ArribaTranslator,
    CiceroTranslator,
    CIVICTranslator,
    EnFusionTranslator,
    FusionCatcherTranslator,
    GenieTranslator,
    JaffaTranslator,
    StarFusionTranslator,
    Translator,
)

_logger = logging.getLogger(__name__)


class FusionCallerHarvester(ABC):
    """ABC for fusion caller harvesters"""

    fusion_caller: FusionCaller
    column_rename: dict
    delimeter: str
    translator: Translator
    coordinate_type: CoordinateType

    def __init__(self, translator: Translator, assembly: Assembly) -> None:
        """Initialize FusionCallerHarvester

        :param translator: A Translator object
        :param assembly: The assembly that the coordinates are described on
        """
        self.translator = translator
        self.assembly = assembly

    def _get_records(self, fusions_file: TextIO) -> csv.DictReader:
        """Read in all records from a fusions file

        :param fusion_file: The open fusions file
        :return A csv.DictReader object containing the detected fusions
        """
        return csv.DictReader(fusions_file, delimiter=self.delimeter)

    @staticmethod
    def _count_dropped_fusions(
        initial_fusions: list,
        translated_fusions: list[AssayedFusion | CategoricalFusion],
    ) -> None:
        """Count the number of fusions that were dropped during translation

        :param initial_fusions: A list containing the original fusions from the provided source
        :param translated_fusions: A list of translated fusions (assayed or categorical)
        """
        diff = len(initial_fusions) - len(translated_fusions)
        if diff:
            msg = f"{diff} fusion(s) were dropped during translation"
            _logger.warning(msg)

    async def load_records(
        self,
        fusion_path: Path,
    ) -> list[FusionCaller]:
        """Convert rows of fusion caller output to Pydantic classes

        :param fusion_path: The path to the fusions file
        :raise ValueError: if the file does not exist at the specified path
        :return: A list of translated fusions, represented as Pydantic objects
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

        translated_fusions = []
        for fusion in fusions_list:
            try:
                translated_fusion = await self.translator.translate(
                    fusion, self.coordinate_type, self.assembly
                )
            except ValueError as error:
                _logger.exception(error)
            else:
                translated_fusions.append(translated_fusion)
        self._count_dropped_fusions(fusions_list, translated_fusions)

        return translated_fusions


class JAFFAHarvester(FusionCallerHarvester):
    """Class for harvesting JAFFA data"""

    column_rename: ClassVar[dict[str, str]] = {
        "fusion genes": "fusion_genes",
        "spanning reads": "spanning_reads",
        "spanning pairs": "spanning_pairs",
    }
    delimeter = ","
    fusion_caller = JAFFA
    coordinate_type = CoordinateType.RESIDUE

    def __init__(self, fusor: FUSOR, assembly: Assembly):
        """Initialze JAFFAHarvester
        :param fusor: A FUSOR object
        :param assembly: The assembly
        """
        super().__init__(translator=JaffaTranslator(fusor), assembly=assembly)


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
    coordinate_type = CoordinateType.RESIDUE

    def __init__(self, fusor: FUSOR, assembly: Assembly):
        """Initialze StarFusionHarvester
        :param fusor: A FUSOR object
        :param assembly: The assembly
        """
        super().__init__(translator=StarFusionTranslator(fusor), assembly=assembly)


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
    coordinate_type = CoordinateType.RESIDUE

    def __init__(self, fusor: FUSOR, assembly: Assembly):
        """Initialze FusionCatcher Harvester
        :param fusor: A FUSOR object
        :param assembly: The assembly
        """
        super().__init__(translator=FusionCatcherTranslator(fusor), assembly=assembly)


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
    coordinate_type = CoordinateType.RESIDUE

    def __init__(self, fusor: FUSOR, assembly: Assembly):
        """Initialze ArribaHarvester
        :param fusor: A FUSOR object
        :param assembly: The assembly
        """
        super().__init__(translator=ArribaTranslator(fusor), assembly=assembly)


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
    coordinate_type = CoordinateType.RESIDUE

    def __init__(self, fusor: FUSOR, assembly: Assembly):
        """Initialze CiceroHarvester
        :param fusor: A FUSOR object
        :param assembly: The assembly
        """
        super().__init__(translator=CiceroTranslator(fusor), assembly=assembly)


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
    coordinate_type = CoordinateType.RESIDUE

    def __init__(self, fusor: FUSOR, assembly: Assembly):
        """Initialze EnFusionHarvester
        :param fusor: A FUSOR object
        :param assembly: The assembly
        """
        super().__init__(translator=EnFusionTranslator(fusor), assembly=assembly)

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
    coordinate_type = CoordinateType.RESIDUE

    def __init__(self, fusor: FUSOR, assembly: Assembly):
        """Initialze GenieHarvester
        :param fusor: A FUSOR object
        :param assembly: The assembly
        """
        super().__init__(translator=GenieTranslator(fusor), assembly=assembly)


class CIVICHarvester(FusionCallerHarvester):
    """Class for harvesting CIViC Fusion objects"""

    def __init__(
        self,
        fusor: FUSOR,
        update_cache: bool = False,
        update_from_remote: bool = True,
        local_cache_path: str = civic.LOCAL_CACHE_PATH,
    ) -> None:
        """Initialize CivicHarvester class.

        :param fusor: A FUSOR class instance
        :param update_cache: ``True`` if civicpy cache should be updated. Note
            this will take several minutes. ``False`` if to use local cache.
        :param update_from_remote: If set to ``True``, civicpy.update_cache will first
            download the remote cache designated by REMOTE_CACHE_URL, store it
            to LOCAL_CACHE_PATH, and then load the downloaded cache into memory.
            This parameter defaults to ``True``.
        :param local_cache_path: A filepath destination for the retrieved remote
            cache. This parameter defaults to LOCAL_CACHE_PATH from civicpy.
        """
        # super().__init__(translator, Assembly.GRCH37)
        if update_cache:
            civic.update_cache(from_remote_cache=update_from_remote)

        civic.load_cache(local_cache_path=local_cache_path, on_stale="ignore")
        self.translator = CIVICTranslator(fusor=fusor)
        self.fusions_list = None

    async def load_records(self) -> list[CIVIC]:
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

        translated_fusions = []
        for fusion in processed_fusions:
            if (
                "?" in fusion.vicc_compliant_name
            ):  # Making suggestion to CIViC to fix syntax (MP: 5474)
                continue
            cat_fusion = await self.translator.translate(civic=fusion)
            translated_fusions.append(cat_fusion)
        self._count_dropped_fusions(processed_fusions, translated_fusions)

        return translated_fusions
