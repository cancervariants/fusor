"""Harvester methods for output from different fusion callers"""

import csv
import json
import logging
from abc import ABC
from collections.abc import AsyncGenerator
from itertools import dropwhile
from pathlib import Path
from typing import Any, ClassVar, Generic, Literal, TextIO, TypeVar

from civicpy import civic
from cool_seq_tool.schemas import Assembly, CoordinateType
from pydantic.dataclasses import dataclass
from wags_tails import MoaData

from fusor.config import config
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
from fusor.models import (
    AbstractTranscriptStructuralVariant,
    AssayedFusion,
    CategoricalFusion,
    InternalTandemDuplication,
)
from fusor.translator import (
    ArribaTranslator,
    CiceroTranslator,
    CIVICTranslator,
    EnFusionTranslator,
    FusionCatcherTranslator,
    GenieTranslator,
    JAFFATranslator,
    MOATranslator,
    STARFusionTranslator,
    Translator,
)

_logger = logging.getLogger(__name__)

T = TypeVar("T", bound=Translator)


@dataclass
class FusionCallerRecord:
    """Records a single entry in a fusion caller data table.

    :param source: The untranslated fusion caller row
    :param annotated: The annotated fusion data
    :param annotation_errors: Captures any errors that occurred during fusion translation
    """

    source: dict[str, Any]
    annotated: AbstractTranscriptStructuralVariant | None
    annotation_error: str | None


class FusionCallerHarvester(ABC, Generic[T]):
    """ABC for fusion caller harvesters"""

    fusion_caller: type[FusionCaller]
    column_rename: ClassVar[dict[str, str]]
    delimiter: str
    translator_class: type[T]
    coordinate_type: CoordinateType

    def __init__(self, fusor: FUSOR, assembly: Assembly) -> None:
        """Initialize FusionCallerHarvester

        :param fusor: A FUSOR object
        :param assembly: The assembly that the coordinates are described on
        """
        self.translator = self.translator_class(fusor)
        self.assembly = assembly

    def _get_records(self, fusions_file: TextIO) -> csv.DictReader:
        """Read in all records from a fusions file

        :param fusion_file: The open fusions file
        :return A csv.DictReader object containing the detected fusions
        """
        return csv.DictReader(fusions_file, delimiter=self.delimiter)

    @staticmethod
    def _count_dropped_fusions(
        initial_fusions: list,
        translated_fusions: list,
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
    ) -> list[AssayedFusion | InternalTandemDuplication]:
        """Convert rows of fusion caller output to AssayedFusion and InternalTandemDuplication objects

        :param fusion_path: The path to the fusions file
        :raise ValueError: if the file does not exist at the specified path
        :return: A list of translated fusions, represented as AssayedFusion or InternalTandemDuplcation objects
        """
        records = await self.load_record_table(fusion_path)
        fusions = [
            record.annotated
            for record in records
            if isinstance(record.annotated, AssayedFusion | InternalTandemDuplication)
        ]
        self._count_dropped_fusions(records, fusions)

        return fusions

    async def load_record_table(
        self,
        fusion_path: Path,
    ) -> list[FusionCallerRecord]:
        """Convert rows of fusion caller output to FusionCallerRecord objects.

        Each entry in the returned list includes the raw fusion caller output row, as well as the
        annotated fusion model. Any errors encountered when annotating fusions are captured and returned per fusion.

        :param fusion_path: The path to the fusions file
        :raise ValueError: if the file does not exist at the specified path
        :return: A list of fusion caller records
        """
        records: list[FusionCallerRecord] = []
        async for raw_row, fusion in self.harvest_records(fusion_path):
            error: str | None = None
            translated_fusion: AbstractTranscriptStructuralVariant | None = None
            try:
                translated_fusion = await self.translate(fusion)
            except Exception as exc:
                _logger.exception(
                    "Error encountered while loading records from %s", fusion_path
                )
                error = str(exc)

            records.append(
                FusionCallerRecord(
                    source=raw_row,
                    annotated=None
                    if isinstance(translated_fusion, str)
                    else translated_fusion,
                    annotation_error=error,
                )
            )

        return records

    async def translate(
        self, fusion: FusionCaller
    ) -> AbstractTranscriptStructuralVariant:
        """Call the translator for this fusion.

        :param fusion: The fusion call entry.
        :return: The translated fusion call. Usually an AssayedFusion.
        """
        return await self.translator.translate(
            fusion, self.coordinate_type, self.assembly
        )

    async def harvest_records(
        self,
        fusion_path: Path,
    ) -> AsyncGenerator[tuple[dict, FusionCaller], None]:
        """Convert rows of fusion caller output to FusionCallerRecord objects.

        Each entry in the returned list includes the raw fusion caller output row, as well as the
        annotated fusion model. Any errors encountered when annotating fusions are captured and returned per fusion.

        :param fusion_path: The path to the fusions file
        :raise ValueError: if the file does not exist at the specified path
        :return: A generator of raw fusion rows, and fusion call outputs.
        """
        if not fusion_path.exists():
            msg = f"{fusion_path!s} does not exist"
            raise ValueError(msg)

        fields_to_keep = self.fusion_caller.__annotations__
        reader = self._get_records(fusion_path.open())

        for row in reader:
            raw_row = {}
            filtered_row = {}
            for key, value in row.items():
                renamed_key = self.column_rename.get(key, key)
                raw_row[renamed_key] = value
                if renamed_key in fields_to_keep:
                    filtered_row[renamed_key] = value

            fusion = self.fusion_caller(**filtered_row)

            yield (raw_row, fusion)


class JAFFAHarvester(FusionCallerHarvester):
    """Class for harvesting JAFFA data"""

    column_rename: ClassVar[dict[str, str]] = {
        "fusion genes": "fusion_genes",
        "spanning reads": "spanning_reads",
        "spanning pairs": "spanning_pairs",
    }
    delimiter = ","
    fusion_caller = JAFFA
    translator_class = JAFFATranslator
    coordinate_type = CoordinateType.RESIDUE


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
    delimiter = "\t"
    fusion_caller = STARFusion
    translator_class = STARFusionTranslator
    coordinate_type = CoordinateType.RESIDUE


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
    delimiter = "\t"
    fusion_caller = FusionCatcher
    translator_class = FusionCatcherTranslator
    coordinate_type = CoordinateType.RESIDUE


class ArribaHarvester(FusionCallerHarvester):
    """Class for harvesting Arriba data"""

    column_rename: ClassVar[dict[str, str]] = {
        "#gene1": "gene1",
        "strand1(gene/fusion)": "strand1",
        "strand2(gene/fusion)": "strand2",
        "type": "event_type",
        "reading_frame": "rf",
    }
    delimiter = "\t"
    fusion_caller = Arriba
    translator_class = ArribaTranslator
    coordinate_type = CoordinateType.RESIDUE


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
    delimiter = "\t"
    fusion_caller = Cicero
    translator_class = CiceroTranslator
    coordinate_type = CoordinateType.RESIDUE


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
    delimiter = "\t"
    fusion_caller = EnFusion
    translator_class = EnFusionTranslator
    coordinate_type = CoordinateType.RESIDUE

    def _get_records(self, fusions_file: TextIO) -> csv.DictReader:
        """Read in all records from a fusions file

        :param fusion_file: The open fusions file
        :return A csv.DictReader object containing the detected fusions
        """
        fusion_lines = dropwhile(
            lambda line: not line.startswith("UnorderedFusion"), fusions_file
        )
        return csv.DictReader(fusion_lines, delimiter=self.delimiter)


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
    delimiter = "\t"
    fusion_caller = Genie
    translator_class = GenieTranslator
    coordinate_type = CoordinateType.RESIDUE


class CIVICHarvester(FusionCallerHarvester):
    """Class for harvesting CIViC Fusion objects"""

    translator_class = CIVICTranslator

    def __init__(
        self,
        fusor: FUSOR,
        update_cache: bool = False,
        update_from_remote: bool = True,
        local_cache_path: str = civic.LOCAL_CACHE_PATH,
        include_status: list[Literal["accepted", "submitted", "rejected"]]
        | None = None,
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
        :param include_status: Whether to include accepted, submitted, and/or
            rejected fusion variants from civicpy cache. By default, this is
            set to ``accepted`` if include_status is set to None.
        """
        super().__init__(fusor, Assembly.GRCH37)
        if update_cache:
            civic.update_cache(from_remote_cache=update_from_remote)

        civic.load_cache(local_cache_path=local_cache_path, on_stale="ignore")
        self.translator = CIVICTranslator(fusor=fusor)

        # Set include_status to accepted if it is None
        if not include_status:
            include_status = ["accepted"]
        # Load in fusion variants from CIViC
        self.fusions_list = civic.get_all_fusion_variants(include_status=include_status)

    async def load_records(self) -> list[CategoricalFusion]:
        """Convert CIViC fusions to CategoricalFusion objects

        :return A list of CategoricalFusion objects
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
            try:
                cat_fusion = await self.translator.translate(civic=fusion)
            except ValueError as e:
                err_msg = f"Cannot translate fusion: {fusion.vicc_compliant_name} due to the following reason: {e!s}"
                _logger.exception(err_msg)
                continue
            if cat_fusion:
                translated_fusions.append(cat_fusion)
        self._count_dropped_fusions(processed_fusions, translated_fusions)

        return translated_fusions


class MOAHarvester(FusionCallerHarvester):
    """Class for harvesting Molecular Oncology Almanac (MOA) fusion data"""

    translator_class: MOATranslator

    def __init__(
        self,
        fusor: FUSOR,
        cache_dir: Path | None = None,
        force_refresh: bool = False,
        use_local: bool = False,
    ) -> None:
        """Initialize MOAHarvester class

        :param fusor: A FUSOR object
        :param cache_dir: The path to the store the cached MOA assertions.
            This by defualt is set to None, and the MOA assertions are
            stored in the `FUSOR_DATA_DIR` directory.
        :paran force_refresh: A boolean indicating if the MOA assertions
            file should be regenerated. By default, this is set to ``False``.
        :param use_local: A boolean indicating if the latest local available
            file should be use. By default, this is set to ``False``.
        """
        self.translator = MOATranslator(fusor)
        if not cache_dir:
            cache_dir = config.data_root
        cache_dir.mkdir(parents=True, exist_ok=True)
        moa_downloader = MoaData(data_dir=cache_dir)
        moa_file = moa_downloader.get_latest(
            force_refresh=force_refresh, from_local=use_local
        )[0]
        with moa_file.open("rb") as f:
            moa_assertions = json.load(f)
            self.assertions = moa_assertions["content"]

    def load_records(self) -> list[CategoricalFusion]:
        """Convert MOA records to CategoricalFusion objects

        :return A list of CategoricalFusion objects
        """
        # Filter assertion dicts to only include fusion events
        moa_fusions = [
            assertion
            for assertion in self.assertions
            if any(
                ext.get("name") == "rearrangement_type" and ext.get("value") == "Fusion"
                for biomarker in assertion.get("proposition", {}).get("biomarkers", [])
                for ext in biomarker.get("extensions", [])
            )
        ]
        translated_fusions = []

        for fusion in moa_fusions:
            moa_fusion = self.translator.translate(fusion)
            if moa_fusion:
                translated_fusions.append(moa_fusion)
        self._count_dropped_fusions(moa_fusions, translated_fusions)

        return translated_fusions
