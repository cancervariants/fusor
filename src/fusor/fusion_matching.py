"""Module for matching assayed fusions against categorical fusions"""

import pickle
from pathlib import Path

from fusor.models import (
    AssayedFusion,
    CategoricalFusion,
    FusionSet,
    GeneElement,
    LinkerElement,
    MultiplePossibleGenesElement,
    TranscriptSegmentElement,
    UnknownGeneElement,
)


class FusionMatcher:
    """Class for matching assayed fusions against categorical fusions"""

    def __init__(
        self,
        cache_dir: Path,
        fusion_set: FusionSet,
        cache_files: list[str] | None = None,
    ) -> None:
        """Initialize FusionMatcher class and comparator categorical fusion objects

        :param cache_dir: The directory containing the cached categorical fusions
            files. If cached files do not exist in the directory, a cached file at
            the provided location will be generated for each source.
        :param fusion_set: A FusionSet object
        :param cache_files: A list of cache file names containing CategoricalFusion
            objects to load, or None. By default this is set to None.
        """
        self.cache_dir = cache_dir
        self.assayed_fusions = fusion_set.assayedFusions
        self.categorical_fusions = fusion_set.categoricalFusions
        if cache_files:
            self.cache_files = cache_files

    async def _load_categorical_fusions(self) -> list[CategoricalFusion]:
        """Load in cache of CategoricalFusion objects

        :return A list of Categorical fusions
        """
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
    ) -> list[CategoricalFusion]:
        """Filter CategoricalFusion list to ensure fusion matching is run on relevant list

        :param assayed_fusion: The AssayedFusion object that is being queried
        :param categorical_fusions: A list of CategoricalFusion objects
        :return A list of filtered categorical fusion objects, or an empty list if no
            filtered categorical fusions are generated
        """
        return [
            categorical_fusion
            for categorical_fusion in categorical_fusions
            if self._match_fusion_partners(assayed_fusion, categorical_fusion)
        ]

    def _match_fusion_structure(
        self,
        assayed_element: TranscriptSegmentElement | UnknownGeneElement | GeneElement,
        categorical_element: TranscriptSegmentElement
        | MultiplePossibleGenesElement
        | GeneElement,
        is_five_prime_partner: bool,
    ) -> int | None:
        """Compare fusion partner information for assayed and categorical fusions. A
            maximum of 5 fields are compared: the gene symbol, transcript accession,
            exon number, exon offset, and genomic breakpoint. A match score of 5 is
            returned if all these fields are equivalent. A match score of 0 is returned
            if one of these fields differs. `None` is returned when a gene partner
            is ambiguous and a clear comparison cannot be made.

        :param assayed_element: The assayed fusion transcript or unknown gene element
            or gene element
        :param categorical_element: The categorical fusion transcript or mulitple
            possible genes element
        :param is_five_prime_partner: If the 5' fusion partner is being compared
        :return A score indiciating the degree of match or None if the partner is
            ? or v
        """
        # Set default match score
        match_score = 0

        # If the assayed partner is unknown or the categorical partner is a multiple
        # possible gene element, return None as no precise information
        # regarding the compared elements is known
        if isinstance(assayed_element, UnknownGeneElement) or isinstance(
            categorical_element, MultiplePossibleGenesElement
        ):
            return None

        # Compare gene partners first
        if assayed_element.gene == categorical_element.gene:
            match_score += 1
        else:
            return 0

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
                return 0

            start_or_end = "End" if is_five_prime_partner else "Start"
            fields_to_compare = [
                f"exon{start_or_end}",
                f"exon{start_or_end}Offset",
                f"elementGenomic{start_or_end}",
            ]

            # Determine if exon number, offset, and genomic breakpoint match
            for field in fields_to_compare:
                if getattr(assayed_element, field) == getattr(
                    categorical_element, field
                ):
                    match_score += 1
                else:
                    return 0

        return match_score

    def _compare_fusion(
        self, assayed_fusion: AssayedFusion, categorical_fusion: CategoricalFusion
    ) -> bool | tuple[bool, int]:
        """Compare assayed and categorical fusions

        :param assayed_fusion: AssayedFusion object
        :param categorical_fusion: CategoricalFusion object
        :return A boolean or a tuple containing a boolean and match score. The boolean
            indicates if the AssayedFusion and CategoricalFusion objects match, and the
            match score describes the quality of the match. A higher match score
            indicates a higher quality match.
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
                return False  # Linker Sequences are not equivalent, no match

        # Compare other structural elements
        match_data_5prime = self._match_fusion_structure(
            assayed_transcript_segments[0], categorical_transcript_segments[0], True
        )
        if match_data_5prime == 0:
            return False
        match_data_3prime = self._match_fusion_structure(
            assayed_transcript_segments[1], categorical_transcript_segments[1], False
        )
        if match_data_3prime == 0:
            return False

        # Update match scores in event partner is ? or v
        match_data_5prime = match_data_5prime if match_data_5prime else 0
        match_data_3prime = match_data_3prime if match_data_3prime else 0
        match_score = match_score + match_data_5prime + match_data_3prime
        return True, match_score if match_score > 0 else False

    async def match_fusion(
        self,
    ) -> list[list[tuple[CategoricalFusion, int]]]:
        """Return best matching fusion

        :return A list of list of tuples containing matching categorical fusion objects
            and their associated match score or None, for each examined AssayedFusion
            object. This method iterates through all supplied AssayedFusion objects to
            find corresponding matches. The match score represents how many attributes
            are equivalent between an AssayedFusion and CategoricalFusion. The
            attributes that are compared include the gene partner, transcript accession,
            exon number, exon offset, and genomic breakpoint. A higher match score
            indicates that more fields were equivalent. Matches are returned for each
            queried AssayedFusion in descending order, with the highest quality match
            reported first.
        """
        matched_fusions = []
        for assayed_fusion in self.assayed_fusions:
            matching_output = []
            categorical_fusions = self._filter_categorical_fusions(
                assayed_fusion,
                self.categorical_fusions
                if self.categorical_fusions
                else await self._load_categorical_fusions(),
            )
            if (
                not categorical_fusions
            ):  # Return empty list if no filtered fusions are generated
                return []

            for categorical_fusion in categorical_fusions:
                match_information = self._compare_fusion(
                    assayed_fusion, categorical_fusion
                )
                if match_information:
                    matching_output.append((categorical_fusion, match_information[1]))
            matched_fusions.append(
                sorted(matching_output, key=lambda x: x[1], reverse=True)
            )

        return matched_fusions
