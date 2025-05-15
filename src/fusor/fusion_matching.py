"""Module for matching assayed fusions against categorical fusions"""

import pickle
from enum import Enum
from pathlib import Path

from civicpy import civic

from fusor.harvester import CIVICHarvester
from fusor.models import (
    AssayedFusion,
    CategoricalFusion,
    GeneElement,
    LinkerElement,
    MultiplePossibleGenesElement,
    TranscriptSegmentElement,
    UnknownGeneElement,
)
from fusor.translator import Translator


class FusionSources(str, Enum):
    """Define CategoricalFusion sources"""

    CIVIC = "CIViC"


class FusionMatcher:
    """Class for matching assayed fusions against categorical fusions"""

    def __init__(
        self,
        translator: Translator,
    ) -> None:
        """Initialize CategoricalFusion source

        :param translator: A translator object
        """
        self.translator = translator
        self.categorical_fusions: list[CategoricalFusion] = []

    def load_categorical_fusions(self, output_dir: Path) -> list[CategoricalFusion]:
        """Load in cache of CategoricalFusion objects

        :param output_dir: A path to the cache
        :return A list of Categorical fusions
        """
        with output_dir.open("rb") as categorical_fusions:
            return pickle.load(categorical_fusions)  # noqa: S301
        self.categorical_fusions = categorical_fusions
        return self.categorical_fusions

    def _extract_fusion_partners(
        self,
        fusion_elements: list[
            UnknownGeneElement
            | MultiplePossibleGenesElement
            | TranscriptSegmentElement
            | LinkerElement
            | GeneElement
        ],
    ) -> list[str, str]:
        """Extract gene symbols for a fusion event to allow for filtering

        :param fusion_elements: A list of possible fusion elements
        :return The two gene symbols involved in the fusion, or ?/v if one partner is not known/provided
        """
        gene_symbols = []
        for element in fusion_elements:
            if isinstance(element, GeneElement | TranscriptSegmentElement):
                gene_symbols.append(element.gene.name)
            elif isinstance(element, UnknownGeneElement):
                gene_symbols.append("?")
            elif isinstance(element, MultiplePossibleGenesElement):
                gene_symbols.append("v")
        return gene_symbols

    def _match_fusion_partners(
        self, assayed_fusion: AssayedFusion, categorical_fusion: CategoricalFusion
    ) -> bool:
        """Determine if assayed fusion and categorical fusion have the same partners

        :param assayed_fusion_gene_symbols: AssayedFusion object
        :param categorical_fusion: CategoricalFusion object
        :return ``True`` if the symbols for the fusion match match, ``False`` if not
        """
        assayed_fusion_gene_symbols = self._extract_fusion_partners(
            assayed_fusion.structure
        )
        categorical_fusion_gene_symbols = self._extract_fusion_partners(
            categorical_fusion.structure
        )
        return (
            assayed_fusion_gene_symbols == categorical_fusion_gene_symbols
            or (
                categorical_fusion_gene_symbols[0] == "v"
                and assayed_fusion_gene_symbols[1] == categorical_fusion_gene_symbols[1]
            )
            or (
                assayed_fusion_gene_symbols[0] == categorical_fusion_gene_symbols[0]
                and categorical_fusion_gene_symbols[1] == "v"
            )
        )

    def _filter_categorical_fusions(
        self,
        assayed_fusion: AssayedFusion,
        categorical_fusions: list[CategoricalFusion],
    ) -> list[CategoricalFusion] | None:
        """Filter CategoricalFusion list to ensure fusion matching is run on relevant list

        :param assayed_fusion: The AssayedFusion object that is being queried
        :param categorical_fusions: A list of CategoricalFusion objects
        :return A list of filtered categorical fusion objects, or None if the list is empty
        """
        fusions_list = [
            categorical_fusion
            for categorical_fusion in categorical_fusions
            if self._match_fusion_partners(assayed_fusion, categorical_fusion)
        ]
        return fusions_list if fusions_list else None

    def _compare_structure(
        self,
        assayed_element: TranscriptSegmentElement | UnknownGeneElement | GeneElement,
        categorical_element: TranscriptSegmentElement
        | MultiplePossibleGenesElement
        | GeneElement,
        is_five_prime_partner: bool,
    ) -> tuple[bool | str, int]:
        """Compare transcript segments for an assayed and categorical fusions
        :param assayed_element: The assayed fusion transcript or unknown gene element or gene element
        :param categorical_element: The categorical fusion transcript or mulitple possible genes element
        :param is_five_prime_partner: If the 5' fusion partner is being compared
        :return A boolean or string indicating if a match is found and a score indiciating the degree of match
        """
        # Set default match score
        match_score = 0

        # If the assayed partner is unknown or the categorical partner is a multiple possible gene element, return match score of 0 as no precise information
        # regarding the compared elements is known
        if isinstance(assayed_element, UnknownGeneElement) or isinstance(
            categorical_element, MultiplePossibleGenesElement
        ):
            return "NA", 0

        # Compare gene partners first
        if assayed_element.gene == categorical_element.gene:
            match_score += 1
        else:
            return False, 0

        # Then compare transcript partners if transcript data exists
        if isinstance(assayed_element, TranscriptSegmentElement) and isinstance(
            categorical_element, TranscriptSegmentElement
        ):
            if (
                assayed_element.transcript
                and categorical_element.transcript
                and assayed_element.transcript == categorical_element.transcript
            ):
                match_score += 1
            else:
                return False, 0

            start_or_end = "End" if is_five_prime_partner else "Start"
            fields_to_compare = [
                f"exon{start_or_end}",
                f"exon{start_or_end}Offset",
                f"elementGenomic{start_or_end}",
            ]

            for field in fields_to_compare:
                if getattr(assayed_element, field) == getattr(
                    categorical_element, field
                ):
                    match_score += 1
                else:
                    return False, 0

        return True, match_score

    def _compare_fusion(
        self, assayed_fusion: AssayedFusion, categorical_fusion: CategoricalFusion
    ) -> bool | tuple[bool, int]:
        """Compare assayed and categorical fusions

        :param assayed_fusion: AssayedFusion object
        :param categorical_fusion: CategoricalFusion object
        :return A boolean or a tuple containing a boolean and match score
        """
        assayed_transcript_segments = assayed_fusion.structure
        categorical_transcript_segments = categorical_fusion.structure
        match_score = 0

        # Check for linker elements first
        if (
            len(assayed_transcript_segments)
            == len(categorical_transcript_segments)
            == 3
        ):
            if assayed_transcript_segments[1] == categorical_transcript_segments[1]:
                match_score += 1
                assayed_transcript_segments.pop(1)
                categorical_transcript_segments.pop(1)
            else:
                return False

        # Compare other structural elements
        match_data_5prime = self._compare_structure(
            assayed_transcript_segments[0], categorical_transcript_segments[0], True
        )
        if not match_data_5prime[0]:
            return False
        match_data_3prime = self._compare_structure(
            assayed_transcript_segments[1], categorical_transcript_segments[1], False
        )
        if not match_data_3prime[0]:
            return False
        return True, match_score + match_data_5prime[1] + match_data_3prime[1]

    async def match_fusion(
        self,
        assayed_fusion: AssayedFusion,
        sources: tuple[FusionSources],
        cache_dir: Path,
    ) -> list[tuple[CategoricalFusion, int]] | None:
        """Return best matching fusion

        :param assayed_fusion: The assayed fusion object
        :param sources: A list of FusionSources
        :param cache_dir: The location where the cached categorical fusion objects are stored
        :return A list of tuples containing matching categorical fusion objects and their associated match score or None
        """
        matched_fusions = []
        categorical_fusion_classes = {FusionSources.CIVIC: CIVICCategoricalFusions}
        sources = {k: v for k, v in categorical_fusion_classes.items() if k in sources}
        for source in sources:
            source_class = sources[source]
            source_class = source_class(self.translator)
            cache_path = cache_dir / source_class.cache_name
            if not cache_path.exists():
                await source_class.save_categorical_fusions(output_dir=cache_path)
            self.categorical_fusions.extend(
                self.load_categorical_fusions(output_dir=cache_path)
            )

        categorical_fusions = self._filter_categorical_fusions(
            assayed_fusion, self.categorical_fusions
        )
        if (
            not categorical_fusions
        ):  # Return none if no applicable cateogorical fusion exist
            return None

        for categorical_fusion in categorical_fusions:
            match_information = self._compare_fusion(assayed_fusion, categorical_fusion)
            if match_information:
                matched_fusions.append((categorical_fusion, match_information[1]))

        self.categorical_fusions.clear()  # Reset categorical fusions list
        return (
            sorted(matched_fusions, key=lambda x: x[1], reverse=True)
            if matched_fusions
            else None
        )


class CIVICCategoricalFusions(FusionMatcher):
    """Class for loading and generating cache of CIViC Categorical Fusions"""

    def __init__(
        self, translator: Translator, cache_name: str = "civic_translated_fusions.pkl"
    ) -> None:
        """Initialize CIVICCategoricalFusions class

        :param translator: A translator object
        :param cache_name: The file name to save the cached CategoricalFusion objects
        """
        # Load in all accepted fusions variants
        super().__init__(translator)
        variants = civic.get_all_fusion_variants(include_status="accepted")
        harvester = CIVICHarvester()
        harvester.fusions_list = variants
        self.fusions_list = harvester.load_records()
        self.cache_name = cache_name

    async def save_categorical_fusions(self, output_dir: Path) -> None:
        """Load fusion variants from CIViC and convert to CategoricalFusion objects

        :param output_dir: The location where the cache will be stored
        """
        civic_fusions = []
        for fusion in self.fusions_list:
            if "?" in fusion.vicc_compliant_name:
                continue
            cex = await self.translator.from_civic(civic=fusion)
            civic_fusions.append(cex)

        output_dir.parent.mkdir(parents=True, exist_ok=True)
        with output_dir.open("wb") as f:
            pickle.dump(civic_fusions, f)
