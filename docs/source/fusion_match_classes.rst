.. fusion_match_classes:

Fusion Match Classes
====================

An important step in interpreting gene fusion events involves matching detected gene
fusions against genomic knowledgebases to extract relevant clinical evidence. FUSOR's
``FusionMatcher`` class performs this operation, returning and classifying different
kinds of matches along with their associated clinical evidence statements.

The list of potential match types is described below.

Exact Match
-----------

``MatchType.EXACT``

This occurs when the gene partner, transcript accession, exon number, exon offset,
and genomic breakpoint for both the 5' and 3' partners are equivalent between the
queried fusion and comparator fusion. This also occurs when the queried and
comparator linker sequences are equivalent, if present.

Shared Genes, 5' Exact
----------------------

``MatchType.SHARED_GENES_FIVE_PRIME_EXACT``

This occurs when both the 5' and 3' gene partners of the queried fusion match those
of the comparator fusion. Additionally, the transcript, exon number, exon offset,
and genomic breakpoint for the 5' partner are equivalent. However, at least one of
these attributes is not equivalent for the 3' partner.

Shared Genes, 3' Exact
----------------------

``MatchType.SHARED_GENES_THREE_PRIME_EXACT``

This occurs when both the 5' and 3' gene partners of the queried fusion match those
of the comparator fusion. Additionally, the transcript, exon number, exon offset,
and genomic breakpoint for the 3' partner are equivalent. However, at least one of
these attributes is not equivalent for the 5' partner.

Shared Genes
------------

``MatchType.SHARED_GENES``

This occurs when only the 5' and 3' gene partners of the queried fusion match those
of the comparator fusion. The other transcript and exon data for the partners is
either not equivalent or not present.

5' Exact
--------

``MatchType.FIVE_PRIME_EXACT``

This occurs when the 5' partner of the queried fusion and comparator fusion have
equivalent gene partners, transcript accessions, exon numbers, exon offsets, and
genomic breakpoints. The 3' partner of the queried fusion is unequal to the 3'
partner of the comparator fusion or is not available.

3' Exact
--------

``MatchType.THREE_PRIME_EXACT``

This occurs when the 3' partner of the queried fusion and comparator fusion have
equivalent gene partners, transcript accessions, exon numbers, exon offsets, and
genomic breakpoints. The 5' partner of the queried fusion is unequal to the 5'
partner of the comparator fusion or is not available.

5' Gene
-------

``MatchType.FIVE_PRIME_GENE``

This occurs when the 5' partner of the queried fusion and comparator fusion have
equivalent gene partners, but have missing or unequal transcript and breakpoint
data. The 3' partner of the queried fusion is unequal to the 3' partner of the
comparator fusion or is not available.

3' Gene
-------

``MatchType.THREE_PRIME_GENE``

This occurs when the 3' partner of the queried fusion and comparator fusion have
equivalent gene partners, but have missing or unequal transcript and breakpoint
data. The 5' partner of the queried fusion is unequal to the 5' partner of the
comparator fusion or is not available.

Fusion Matching Criteria and Example
------------------------------------

A visualization describing the match criteria is shown below. An example assayed to
categorical fusion match is also shown, where the assayed fusion describes joining
of exon 13 of ``EML4`` with exon 20 of ``ALK``, and the categorical fusion describes
the joining of exon 20 of ``EML4`` with exon 20 of ``ALK``.

In this example, both fusions have the same gene partners, there is an exact match
for the 3' partner, but the 5' exon and 5' junction location are unequal, resulting
in a ``SHARED_GENES_THREE_PRIME_EXACT`` match being returned.

.. figure:: ../../fusion_match_classes.png
   :alt: Fusion Match Classes
   :align: center
