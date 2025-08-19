"""Module for matching assayed fusions against categorical fusions"""

import pickle
from enum import Enum, unique
from pathlib import Path

from pydantic import BaseModel

from fusor.config import config
from fusor.models import (
    AssayedFusion,
    CategoricalFusion,
    GeneElement,
    LinkerElement,
    MultiplePossibleGenesElement,
    TranscriptSegmentElement,
    UnknownGeneElement,
)


@unique
class MatchType(str, Enum):
    """Enum for defining different match types"""

    EXACT = "EXACT"
    SHARED_GENES_FIVE_PRIME_EXACT = "SHARED_GENES_FIVE_PRIME_EXACT"
    SHARED_GENES_THREE_PRIME_EXACT = "SHARED_GENES_THREE_PRIME_EXACT"
    SHARED_GENES = "SHARED_GENES"
    FIVE_PRIME_GENE = "FIVE_PRIME_GENE"
    FIVE_PRIME_EXACT = "FIVE_PRIME_EXACT"
    THREE_PRIME_GENE = "THREE_PRIME_GENE"
    THREE_PRIME_EXACT = "THREE_PRIME_EXACT"
    NO_MATCH = "NO_MATCH"

    @property
    def priority(self) -> int:
        """Return numeric priority for sorting, where the lower the score
        indicates the higher quality match
        """
        return {
            MatchType.EXACT: 10,
            MatchType.SHARED_GENES_FIVE_PRIME_EXACT: 20,
            MatchType.SHARED_GENES_THREE_PRIME_EXACT: 21,
            MatchType.SHARED_GENES: 30,
            MatchType.FIVE_PRIME_EXACT: 40,
            MatchType.THREE_PRIME_EXACT: 41,
            MatchType.FIVE_PRIME_GENE: 50,
            MatchType.THREE_PRIME_GENE: 51,
            MatchType.NO_MATCH: 60,
        }[self]


class MatchInformation(BaseModel):
    """Class for reporting matching information"""

    five_prime_gene: bool | None = None
    five_prime_transcript: bool | None = None
    five_prime_exon: bool | None = None
    five_prime_exon_offset: bool | None = None
    five_prime_breakpoint: bool | None = None
    linker: bool | None = None
    three_prime_gene: bool | None = None
    three_prime_transcript: bool | None = None
    three_prime_exon: bool | None = None
    three_prime_exon_offset: bool | None = None
    three_prime_breakpoint: bool | None = None

    def _transcript_match(self, transcript_data: list[bool | None]) -> bool:
        """Determine if transcript data matches

        :param transcript_data: A list containing transcript elements
        :return: ``True`` if match, ``False if not
        """
        return all(transcript_data)

    def determine_match(self) -> MatchType:
        """Determine match type based on fields in MatchInformation class

        :return: A MatchType object
        """
        five_prime = [
            self.five_prime_gene,
            self.five_prime_transcript,
            self.five_prime_exon,
            self.five_prime_exon_offset,
            self.five_prime_breakpoint,
        ]
        three_prime = [
            self.three_prime_gene,
            self.three_prime_transcript,
            self.three_prime_exon,
            self.three_prime_exon_offset,
            self.three_prime_breakpoint,
        ]

        # Define and return match criteria
        if (
            self._transcript_match(five_prime)
            and self._transcript_match(three_prime)
            and self.linker is None
        ):
            return MatchType.EXACT
        if (
            self._transcript_match(five_prime)
            and self._transcript_match(three_prime)
            and self.linker
        ):
            return MatchType.EXACT
        if (
            self._transcript_match(five_prime)
            and self.three_prime_gene
            and not self._transcript_match(three_prime)
        ):
            return MatchType.SHARED_GENES_FIVE_PRIME_EXACT
        if (
            self._transcript_match(three_prime)
            and self.five_prime_gene
            and not self._transcript_match(five_prime)
        ):
            return MatchType.SHARED_GENES_THREE_PRIME_EXACT
        if (
            self.five_prime_gene
            and self.three_prime_gene
            and not self._transcript_match(five_prime)
            and not self._transcript_match(three_prime)
        ):
            return MatchType.SHARED_GENES
        if self._transcript_match(five_prime) and not self._transcript_match(
            three_prime
        ):
            return MatchType.FIVE_PRIME_EXACT
        if self._transcript_match(three_prime) and not self._transcript_match(
            five_prime
        ):
            return MatchType.THREE_PRIME_EXACT
        if (
            self.five_prime_gene
            and not self._transcript_match(five_prime)
            and not self._transcript_match(three_prime)
        ):
            return MatchType.FIVE_PRIME_GENE
        if (
            self.three_prime_gene
            and not self._transcript_match(three_prime)
            and not self._transcript_match(five_prime)
        ):
            return MatchType.THREE_PRIME_GENE
        return MatchType.NO_MATCH


class FusionMatcher:
    """Class for matching assayed fusions against categorical fusions"""

    def __init__(
        self,
        cache_dir: Path | None = None,
        assayed_fusions: list[AssayedFusion] | None = None,
        categorical_fusions: list[CategoricalFusion] | None = None,
        cache_files: list[str] | None = None,
    ) -> None:
        """Initialize FusionMatcher class and comparator categorical fusion objects

        :param cache_dir: The directory containing the cached categorical fusions
            files. If this parameter is not provided, it will be set by default
            to be `FUSOR_DATA_DIR`.
        :param assayed_fusions: A list of AssayedFusion objects
        :param categorical_fusions: A list of CategoricalFusion objects
        :param cache_files: A list of cache file names in ``cache_dir`` containing
            CategoricalFusion objects to load, or None. By default this is set to None.
            It assumes that files contain lists of valid CategoricalFusion objects.
        :raises ValueError: If ``categorical_fusions`` is not provided and either
            ``cache_dir`` or ``cache_files`` is not provided.
        """
        if not categorical_fusions and (not cache_dir or not cache_files):
            msg = "Either a list of CategoricalFusion objects must be provided to `categorical_fusions` or a Path and list of file names must be provided to `cache_dir` and `cache_files`, respectively"
            raise ValueError(msg)
        if not cache_dir:
            cache_dir = config.data_root
        cache_dir.mkdir(parents=True, exist_ok=True)
        self.cache_dir = cache_dir
        self.assayed_fusions = assayed_fusions
        self.cache_files = cache_files

        # Load in CategoricalFusions, prioritizing those directly provided by the user
        # with self.categorical_fusions
        self.categorical_fusions = (
            categorical_fusions
            if categorical_fusions
            else self._load_categorical_fusions()
        )

    def _load_categorical_fusions(self) -> list[CategoricalFusion]:
        """Load in cache of CategoricalFusion objects

        :raises ValueError: If the cache_dir or cache_files variables are None
        :return: A list of Categorical fusions
        """
        if not self.cache_dir or not self.cache_files:
            msg = "`cache_dir` and `cache_files` parameters must be provided"
            raise ValueError(msg)
        categorical_fusions = []
        for file in self.cache_files:
            cached_file = self.cache_dir / file
            if cached_file.is_file():
                with cached_file.open("rb") as f:
                    categorical_fusions.extend(pickle.load(f))  # noqa: S301
        return categorical_fusions

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
        :return: The two gene symbols involved in the fusion, or ?/v if one partner is not known/provided
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

        :param assayed_fusion: AssayedFusion object
        :param categorical_fusion: CategoricalFusion object
        :return: ``True`` if the symbols for the fusion match, ``False`` if not
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
    ) -> list[CategoricalFusion]:
        """Filter CategoricalFusion list to ensure fusion matching is run on fusions
        whose partners match those in the AssayedFusion object

        :param assayed_fusion: The AssayedFusion object that is being queried
        :return: A list of filtered categorical fusion objects, or an empty list if no
            filtered categorical fusions are generated
        """
        return [
            categorical_fusion
            for categorical_fusion in self.categorical_fusions
            if self._match_fusion_partners(assayed_fusion, categorical_fusion)
        ]

    def _match_fusion_structure(
        self,
        assayed_element: TranscriptSegmentElement | UnknownGeneElement | GeneElement,
        categorical_element: TranscriptSegmentElement
        | MultiplePossibleGenesElement
        | GeneElement,
        is_five_prime_partner: bool,
        mi: MatchInformation,
    ) -> None:
        """Compare fusion partner information for assayed and categorical fusions. A
        maximum of 5 fields are compared: the gene symbol, transcript accession,
        exon number, exon offset, and genomic breakpoint. The supplied
        MatchInformation object is updated throughout the function.

        :param assayed_element: The assayed fusion transcript or unknown gene element
            or gene element
        :param categorical_element: The categorical fusion transcript or mulitple
            possible genes element
        :param is_five_prime_partner: If the 5' fusion partner is being compared
        :param mi: A MatchInformation object
        """
        # If the assayed partner is unknown or the categorical partner is a multiple
        # possible gene element, return None as no precise information
        # regarding the compared elements is known
        if isinstance(assayed_element, UnknownGeneElement) or isinstance(
            categorical_element, MultiplePossibleGenesElement
        ):
            return

        # Compare gene partners first
        if assayed_element.gene == categorical_element.gene:
            if is_five_prime_partner:
                mi.five_prime_gene = True
            else:
                mi.three_prime_gene = True

        # Then compare transcript partners if transcript data exists
        if isinstance(assayed_element, TranscriptSegmentElement) and isinstance(
            categorical_element, TranscriptSegmentElement
        ):
            if assayed_element.transcript == categorical_element.transcript:
                if is_five_prime_partner:
                    mi.five_prime_transcript = True
                else:
                    mi.three_prime_transcript = True

            start_or_end = "End" if is_five_prime_partner else "Start"
            fields_to_compare = [
                ("exon", f"exon{start_or_end}"),
                ("exon_offset", f"exon{start_or_end}Offset"),
                ("breakpoint", f"elementGenomic{start_or_end}"),
            ]

            # Determine if exon number, offset, and genomic breakpoint match
            for mi_field, element_field in fields_to_compare:
                if getattr(assayed_element, element_field) == getattr(
                    categorical_element, element_field
                ):
                    if is_five_prime_partner:
                        setattr(mi, f"five_prime_{mi_field}", True)
                    else:
                        setattr(mi, f"three_prime_{mi_field}", True)

    def _compare_fusion(
        self, assayed_fusion: AssayedFusion, categorical_fusion: CategoricalFusion
    ) -> MatchType:
        """Compare assayed and categorical fusions to determine if their attributes
        are equivalent. If one attribute does not match, then we know the fusions
        do not match.

        :param assayed_fusion: AssayedFusion object
        :param categorical_fusion: CategoricalFusion object
        :return: A MatchType object reporting the type of match
        """
        assayed_fusion_structure = assayed_fusion.structure
        categorical_fusion_structure = categorical_fusion.structure
        match_info = MatchInformation()

        # Check for linker elements first
        if (
            len(assayed_fusion_structure) == len(categorical_fusion_structure) == 3
            and assayed_fusion_structure[1] == categorical_fusion_structure[1]
        ):
            match_info.linker = True
            # Remove linker sequences for additional comparison
            assayed_fusion_structure.pop(1)
            categorical_fusion_structure.pop(1)

        # Compare other structural elements
        self._match_fusion_structure(
            assayed_fusion_structure[0],
            categorical_fusion_structure[0],
            True,
            match_info,
        )
        self._match_fusion_structure(
            assayed_fusion_structure[1],
            categorical_fusion_structure[1],
            False,
            match_info,
        )

        # Determine and return match type
        return match_info.determine_match()

    async def match_fusion(
        self,
    ) -> list[list[tuple[CategoricalFusion, MatchType]] | None]:
        """Return best matching fusion

        This method prioritizes using categorical fusion objects that are
        provided in ``self.categorical_fusions`` as opposed those that exist in the
        ``cache_dir`` directory.

        :raises ValueError: If a list of AssayedFusion objects is not provided
        :return: A list of list of tuples containing matching categorical fusion objects
            and their associated match type, for each examined AssayedFusion
            object. This method iterates through all supplied AssayedFusion objects to
            find corresponding matches. The match type represents how many attributes
            are shared between an AssayedFusion and CategoricalFusion. The
            attributes that are compared include the gene partner, transcript accession,
            exon number, exon offset, and genomic breakpoint. Matches are returned
            according to the priority of their match type.
        """
        if not self.assayed_fusions:
            msg = "`assayed_fusions` must be provided a list of AssayedFusion objects before running `match_fusion`"
            raise ValueError(msg)
        matched_fusions = []

        for assayed_fusion in self.assayed_fusions:
            matching_output = []
            filtered_categorical_fusions = self._filter_categorical_fusions(
                assayed_fusion,
            )
            if not filtered_categorical_fusions:  # Add None to matched_fusions
                matched_fusions.append(None)
                continue

            for categorical_fusion in filtered_categorical_fusions:
                match_type = self._compare_fusion(assayed_fusion, categorical_fusion)
                matching_output.append((categorical_fusion, match_type))
            matched_fusions.append(sorted(matching_output, key=lambda x: x[1].priority))

        return matched_fusions
