"""Module for matching assayed fusions against categorical fusions"""

import pickle
from pathlib import Path

from fusor.models import (
    AssayedFusion,
    CategoricalFusion,
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
        cache_dir: Path | None = None,
        assayed_fusions: list[AssayedFusion] | None = None,
        categorical_fusions: list[CategoricalFusion] | None = None,
        cache_files: list[str] | None = None,
    ) -> None:
        """Initialize FusionMatcher class and comparator categorical fusion objects

        :param cache_dir: The directory containing the cached categorical fusions
            files. If cached files do not exist in the directory, a cached file at
            the provided location will be generated for each source.
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
        :return: A score indiciating the degree of match or None if the partner is
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
            if assayed_element.transcript == categorical_element.transcript:
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
    ) -> tuple[bool, int]:
        """Compare assayed and categorical fusions to determine if their attributes
        are equivalent. If one attribute does not match, then we know the fusions
        do not match.

        :param assayed_fusion: AssayedFusion object
        :param categorical_fusion: CategoricalFusion object
        :return: A tuple containing a boolean and match score. The boolean
            indicates if the AssayedFusion and CategoricalFusion objects match, and the
            match score describes the quality of the match. A higher match score
            indicates a higher quality match.
        """
        assayed_fusion_structure = assayed_fusion.structure
        categorical_fusion_structure = categorical_fusion.structure
        match_score = 0

        # Check for linker elements first
        if len(assayed_fusion_structure) == len(categorical_fusion_structure) == 3:
            if assayed_fusion_structure[1] == categorical_fusion_structure[1]:
                match_score += 1
                # Remove linker sequences for additional comparison
                assayed_fusion_structure.pop(1)
                categorical_fusion_structure.pop(1)
            else:
                return False, 0  # Linker Sequences are not equivalent, no match

        # Compare other structural elements
        match_data_5prime = self._match_fusion_structure(
            assayed_fusion_structure[0], categorical_fusion_structure[0], True
        )
        if match_data_5prime == 0:
            return False, 0
        match_data_3prime = self._match_fusion_structure(
            assayed_fusion_structure[1], categorical_fusion_structure[1], False
        )
        if match_data_3prime == 0:
            return False, 0

        # Update match scores in event partner is ? or v
        match_data_5prime = match_data_5prime if match_data_5prime else 0
        match_data_3prime = match_data_3prime if match_data_3prime else 0
        match_score = match_score + match_data_5prime + match_data_3prime
        return (True, match_score) if match_score > 0 else (False, 0)

    async def match_fusion(
        self,
    ) -> list[list[tuple[CategoricalFusion, int]] | None]:
        """Return best matching fusion

        This method prioritizes using categorical fusion objects that are
        provided in ``self.categorical_fusions`` as opposed those that exist in the
        ``cache_dir`` directory.

        :raises ValueError: If a list of AssayedFusion objects is not provided
        :return: A list of list of tuples containing matching categorical fusion objects
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
                is_match, match_score = self._compare_fusion(
                    assayed_fusion, categorical_fusion
                )
                if is_match:
                    matching_output.append((categorical_fusion, match_score))
            matched_fusions.append(
                sorted(matching_output, key=lambda x: x[1], reverse=True)
            )

        return matched_fusions
