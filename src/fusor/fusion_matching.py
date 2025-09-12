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
        }[self]


class PartnerMatch(BaseModel):
    """Class for describing matching fields for a fusion partner"""

    gene: bool = False
    transcript: bool = False
    exon: bool = False
    exon_offset: bool = False
    breakpoint: bool = False

    @property
    def partner_match(self) -> bool:
        """Check that all values are set to True

        :return: ``True`` if all fields are True, ``False`` if not
        """
        return all(self.__dict__.values())


class MatchInformation(BaseModel):
    """Helper for reporting matching information based off of MatchType"""

    five_prime_match_info: PartnerMatch
    linker: bool | None = None
    three_prime_match_info: PartnerMatch

    def determine_match(self) -> MatchType | None:
        """Determine match type based on fields in MatchInformation class

        :return: A MatchType object, or None if no match exists
        """
        five_prime_match = self.five_prime_match_info.partner_match
        three_prime_match = self.three_prime_match_info.partner_match

        # Define and return match criteria
        if (
            five_prime_match
            and three_prime_match
            and self.linker is None  # Consider exact match if linker is not provided
        ):
            return MatchType.EXACT
        if five_prime_match and three_prime_match and self.linker:
            return MatchType.EXACT

        if (
            five_prime_match
            and self.three_prime_match_info.gene
            and not three_prime_match
        ):
            return MatchType.SHARED_GENES_FIVE_PRIME_EXACT

        if (
            three_prime_match
            and self.five_prime_match_info.gene
            and not five_prime_match
        ):
            return MatchType.SHARED_GENES_THREE_PRIME_EXACT

        if (
            self.five_prime_match_info.gene
            and self.three_prime_match_info.gene
            and not five_prime_match
            and not three_prime_match
        ):
            return MatchType.SHARED_GENES

        if five_prime_match and not three_prime_match:
            return MatchType.FIVE_PRIME_EXACT

        if three_prime_match and not five_prime_match:
            return MatchType.THREE_PRIME_EXACT

        if (
            self.five_prime_match_info.gene
            and not five_prime_match
            and not three_prime_match
        ):
            return MatchType.FIVE_PRIME_GENE

        if (
            self.three_prime_match_info.gene
            and not three_prime_match
            and not five_prime_match
        ):
            return MatchType.THREE_PRIME_GENE
        return None


class FusionMatcher:
    """Class for matching assayed fusions against assayed fusions and categorical fusions"""

    def __init__(
        self,
        cache_dir: Path | None = None,
        assayed_fusions: list[AssayedFusion] | None = None,
        comparator_fusions: list[AssayedFusion | CategoricalFusion] | None = None,
        cache_files: list[str] | None = None,
    ) -> None:
        """Initialize FusionMatcher class and comparator categorical fusion objects

        :param cache_dir: The directory containing the cached categorical fusions
            files. If this parameter is not provided, it will be set by default
            to be `FUSOR_DATA_DIR`.
        :param assayed_fusions: A list of AssayedFusion objects
        :param comparator_fusions: A list of AssayedFusion or CategoricalFusion objects
        :param cache_files: A list of cache file names in ``cache_dir`` containing
            AssayedFusion or CategoricalFusion objects to load, or None. By
            default this is set to None. It assumes that files contain lists
            of valid AssayedFusion or CategoricalFusion objects.
        :raises ValueError: If ``comparator_fusions`` is not provided and either
            ``cache_dir`` or ``cache_files`` is not provided.
        """
        if not comparator_fusions and (not cache_dir or not cache_files):
            msg = "Either a list of comparator fusion objects must be provided to `comparator_fusions` or a Path and list of file names must be provided to `cache_dir` and `cache_files`, respectively"
            raise ValueError(msg)
        if not cache_dir:
            cache_dir = config.data_root
        cache_dir.mkdir(parents=True, exist_ok=True)
        self.cache_dir = cache_dir
        self.assayed_fusions = assayed_fusions
        self.cache_files = cache_files

        # Load in comparator fusions prioritizing those directly provided by the user
        # with self.categorical_fusions
        self.comparator_fusions = (
            comparator_fusions
            if comparator_fusions
            else self._load_comparator_fusions()
        )

    def _load_comparator_fusions(self) -> list[AssayedFusion | CategoricalFusion]:
        """Load in cache of AssayedFusion or CategoricalFusion objects

        :raises ValueError: If the cache_dir or cache_files variables are None
        :return: A list of AssayedFusions or CategoricalFusions
        """
        if not self.cache_dir or not self.cache_files:
            msg = "`cache_dir` and `cache_files` parameters must be provided"
            raise ValueError(msg)
        comparator_fusions = []
        for file in self.cache_files:
            cached_file = self.cache_dir / file
            if cached_file.is_file():
                with cached_file.open("rb") as f:
                    comparator_fusions.extend(pickle.load(f))  # noqa: S301
        return comparator_fusions

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
        self,
        assayed_fusion: AssayedFusion,
        comparator_fusion: AssayedFusion | CategoricalFusion,
    ) -> bool:
        """Determine if assayed fusion and categorical fusion have the shared partners

        :param assayed_fusion: AssayedFusion object
        :param comparator_fusion: AssayedFusion or CategoricalFusion object
        :return: ``True`` if at least one symbol is shared, ``False`` if not
        """
        assayed_fusion_gene_symbols = self._extract_fusion_partners(
            assayed_fusion.structure
        )
        comparator_fusion_gene_symbols = self._extract_fusion_partners(
            comparator_fusion.structure
        )
        return bool(
            set(assayed_fusion_gene_symbols) and set(comparator_fusion_gene_symbols)
        )

    def _filter_comparator_fusions(
        self,
        assayed_fusion: AssayedFusion,
    ) -> list[AssayedFusion | CategoricalFusion]:
        """Filter comparator list to ensure fusion matching is run on fusions
        whose partners match those in the AssayedFusion object

        :param assayed_fusion: The AssayedFusion object that is being queried
        :return: A list of filtered AssayedFusion or CategoricalFusion, or an
            empty list if no filtered fusions are generated
        """
        return [
            comparator_fusion
            for comparator_fusion in self.comparator_fusions
            if self._match_fusion_partners(assayed_fusion, comparator_fusion)
        ]

    def _match_fusion_structure(
        self,
        assayed_element: TranscriptSegmentElement | UnknownGeneElement | GeneElement,
        comparator_element: TranscriptSegmentElement
        | MultiplePossibleGenesElement
        | UnknownGeneElement
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
        :param comparator_element: The comparator fusion transcript,
            MultiplePossibleGenesElement, UnknownGeneElement, or GeneElement
        :param is_five_prime_partner: If the 5' fusion partner is being compared
        :param mi: A MatchInformation object
        """
        # If the assayed partner is unknown or the categorical partner is a multiple
        # possible gene element or unknown, return None as no precise information
        # regarding the compared elements is known
        if isinstance(assayed_element, UnknownGeneElement) or isinstance(
            comparator_element, MultiplePossibleGenesElement | UnknownGeneElement
        ):
            return

        # Compare gene partners first
        if assayed_element.gene == comparator_element.gene:
            if is_five_prime_partner:
                mi.five_prime_match_info.gene = True
            else:
                mi.three_prime_match_info.gene = True

        # Then compare transcript partners if transcript data exists
        if isinstance(assayed_element, TranscriptSegmentElement) and isinstance(
            comparator_element, TranscriptSegmentElement
        ):
            if assayed_element.transcript == comparator_element.transcript:
                if is_five_prime_partner:
                    mi.five_prime_match_info.transcript = True
                else:
                    mi.three_prime_match_info.transcript = True

            start_or_end = "End" if is_five_prime_partner else "Start"
            fields_to_compare = [
                ("exon", f"exon{start_or_end}"),
                ("exon_offset", f"exon{start_or_end}Offset"),
                ("breakpoint", f"elementGenomic{start_or_end}"),
            ]

            # Determine if exon number, offset, and genomic breakpoint match
            for mi_field, element_field in fields_to_compare:
                if getattr(assayed_element, element_field) == getattr(
                    comparator_element, element_field
                ):
                    if is_five_prime_partner:
                        setattr(mi.five_prime_match_info, mi_field, True)
                    else:
                        setattr(mi.three_prime_match_info, mi_field, True)

    def _compare_fusion(
        self,
        assayed_fusion: AssayedFusion,
        comparator_fusion: AssayedFusion | CategoricalFusion,
    ) -> MatchType | None:
        """Compare assayed and categorical fusions to determine if their attributes
        are equivalent. If one attribute does not match, then we know the fusions
        do not match.

        :param assayed_fusion: AssayedFusion object
        :param comparator_fusion: AssayedFusion or CategoricalFusion object
        :return: A MatchType object reporting the type of match, or None if no
            match exists
        """
        assayed_fusion_structure = assayed_fusion.structure
        comparator_fusion_structure = comparator_fusion.structure
        match_info = MatchInformation(
            five_prime_match_info=PartnerMatch(), three_prime_match_info=PartnerMatch()
        )

        # Check for linker elements first
        if len(assayed_fusion_structure) == len(comparator_fusion_structure) == 3:
            match_info.linker = (
                assayed_fusion_structure[1] == comparator_fusion_structure[1]
            )
            # Remove linker sequences for additional comparison
            assayed_fusion_structure.pop(1)
            comparator_fusion_structure.pop(1)

        # Compare other structural elements
        self._match_fusion_structure(
            assayed_fusion_structure[0],
            comparator_fusion_structure[0],
            True,
            match_info,
        )
        self._match_fusion_structure(
            assayed_fusion_structure[1],
            comparator_fusion_structure[1],
            False,
            match_info,
        )

        # Determine and return match type
        return match_info.determine_match()

    async def match_fusion(
        self,
    ) -> list[list[tuple[AssayedFusion | CategoricalFusion, MatchType]] | None]:
        """Return best matching fusion

        This method prioritizes using categorical fusion objects that are
        provided in ``self.categorical_fusions`` as opposed those that exist in the
        ``cache_dir`` directory.

        :raises ValueError: If a list of AssayedFusion objects is not provided
        :return: A list of list of tuples containing matching fusion objects
            and their associated match type, for each examined AssayedFusion
            object. This method iterates through all supplied AssayedFusion objects to
            find corresponding matches. The match type represents how many attributes
            are shared between an AssayedFusion and comparator fusion. The
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
            filtered_comparator_fusions = self._filter_comparator_fusions(
                assayed_fusion,
            )
            if not filtered_comparator_fusions:  # Add None to matched_fusions
                matched_fusions.append(None)
                continue

            for comparator_fusion in filtered_comparator_fusions:
                match_type = self._compare_fusion(assayed_fusion, comparator_fusion)
                if match_type:  # Add comparator fusion if there is a match
                    matching_output.append((comparator_fusion, match_type))
            matched_fusions.append(sorted(matching_output, key=lambda x: x[1].priority))

        return matched_fusions
