"""Module for testing FUSOR Translators"""

import polars as pl
import pytest
from cool_seq_tool.schemas import Assembly, CoordinateType

from fusor.fusion_caller_models import (
    JAFFA,
    Arriba,
    Caller,
    Cicero,
    EnFusion,
    FusionCatcher,
    Genie,
    STARFusion,
)
from fusor.models import (
    AnchoredReads,
    AssayedFusion,
    BreakpointCoverage,
    ContigSequence,
    ReadData,
    SpanningReads,
    SplitReads,
    UnknownGeneElement,
)


@pytest.fixture(scope="module")
def fusion_data_example():
    """Create example assayed fusion for TPM3::PDGFRB with exonic breakpoints"""

    def _create_base_fixture(**kwargs):
        params = {
            "type": "AssayedFusion",
            "structure": [
                {
                    "type": "TranscriptSegmentElement",
                    "transcript": "refseq:NM_152263.4",
                    "exonEnd": 8,
                    "exonEndOffset": -66,
                    "gene": {
                        "concept_id": "hgnc:12012",
                        "type": "Gene",
                        "symbol": "TPM3",
                    },
                    "elementGenomicEnd": {
                        "id": "ga4gh:SL.6lXn5i3zqcZUfmtBSieTiVL4Nt2gPGKY",
                        "type": "SequenceLocation",
                        "digest": "6lXn5i3zqcZUfmtBSieTiVL4Nt2gPGKY",
                        "sequenceReference": {
                            "id": "refseq:NC_000001.11",
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        },
                        "start": 154170465,
                    },
                },
                {
                    "type": "TranscriptSegmentElement",
                    "transcript": "refseq:NM_002609.4",
                    "exonStart": 11,
                    "exonStartOffset": 2,
                    "gene": {
                        "concept_id": "hgnc:8804",
                        "type": "Gene",
                        "symbol": "PDGFRB",
                    },
                    "elementGenomicStart": {
                        "id": "ga4gh:SL.Sp1lwuHbRCkWIoe4zzwVKPsS8zK8i0ck",
                        "type": "SequenceLocation",
                        "digest": "Sp1lwuHbRCkWIoe4zzwVKPsS8zK8i0ck",
                        "sequenceReference": {
                            "id": "refseq:NC_000005.10",
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI",
                        },
                        "end": 150126612,
                    },
                },
            ],
            "causativeEvent": {"type": "CausativeEvent", "eventType": "rearrangement"},
            "r_frame_preserved": True,
            "assay": None,
        }
        assayed_fusion = AssayedFusion(**params)
        return assayed_fusion.model_copy(update=kwargs)

    return _create_base_fixture


@pytest.fixture(scope="module")
def fusion_data_example_nonexonic():
    """Create example assayed fusion for TPM3::PDGFRB with non-exonic breakpoints"""

    def _create_base_fixture(**kwargs):
        params = {
            "type": "AssayedFusion",
            "structure": [
                {
                    "type": "TranscriptSegmentElement",
                    "transcript": "refseq:NM_152263.4",
                    "exonEnd": 4,
                    "exonEndOffset": 5,
                    "gene": {
                        "concept_id": "hgnc:12012",
                        "type": "Gene",
                        "symbol": "TPM3",
                    },
                    "elementGenomicEnd": {
                        "id": "ga4gh:SL.O1rVKQA2FTdy_FFWg3qJVSTG_TF_Mkex",
                        "type": "SequenceLocation",
                        "digest": "O1rVKQA2FTdy_FFWg3qJVSTG_TF_Mkex",
                        "sequenceReference": {
                            "id": "refseq:NC_000001.11",
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        },
                        "start": 154173078,
                    },
                },
                {
                    "type": "TranscriptSegmentElement",
                    "transcript": "refseq:NM_002609.4",
                    "exonStart": 11,
                    "exonStartOffset": -559,
                    "gene": {
                        "concept_id": "hgnc:8804",
                        "type": "Gene",
                        "symbol": "PDGFRB",
                    },
                    "elementGenomicStart": {
                        "id": "ga4gh:SL.GtoWMuox4tOyX2I5L9Baobnpgc1pDIVJ",
                        "type": "SequenceLocation",
                        "digest": "GtoWMuox4tOyX2I5L9Baobnpgc1pDIVJ",
                        "sequenceReference": {
                            "id": "refseq:NC_000005.10",
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI",
                        },
                        "end": 150127173,
                    },
                },
            ],
            "causativeEvent": {"type": "CausativeEvent", "eventType": "rearrangement"},
            "r_frame_preserved": True,
            "assay": None,
        }
        assayed_fusion = AssayedFusion(**params)
        return assayed_fusion.model_copy(update=kwargs)

    return _create_base_fixture


def test_gene_element_arriba(translator_instance):
    """Test gene selection for Arriba"""
    genes = "RP1-222H5.1(151985),MIR3672(13973)"
    gene = translator_instance._get_gene_element(genes=genes, caller=Caller.ARRIBA)
    assert gene.gene.symbol == "MIR3672"


def test_valid_fusion_partners(translator_instance):
    """Test that the fusion partners supplied to the translator are different"""
    partners_check = translator_instance._are_fusion_partners_different("BCR", "ABL1")
    assert partners_check

    partners_check = translator_instance._are_fusion_partners_different("BCR", "BCR")
    assert not partners_check


@pytest.mark.asyncio()
async def test_jaffa(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test JAFFA translator"""
    # Test exonic breakpoint
    jaffa = JAFFA(
        fusion_genes="TPM3:PDGFRB",
        chrom1="chr1",
        base1=154170465,
        chrom2="chr5",
        base2=150126612,
        rearrangement=True,
        classification="HighConfidence",
        inframe=True,
        spanning_reads=100,
        spanning_pairs=80,
    )

    jaffa_fusor = await translator_instance.from_jaffa(
        jaffa,
        CoordinateType.INTER_RESIDUE.value,
        Assembly.GRCH38.value,
    )
    fusion_data_example = fusion_data_example(
        readData=ReadData(
            split=SplitReads(splitReads=100), spanning=SpanningReads(spanningReads=80)
        )
    )
    assert jaffa_fusor.structure == fusion_data_example.structure
    assert jaffa_fusor.readData == fusion_data_example.readData

    # Test non-exonic breakpoint
    jaffa.base1 = 154173079
    jaffa.base2 = 150127173

    jaffa_fusor_nonexonic = await translator_instance.from_jaffa(
        jaffa,
        CoordinateType.RESIDUE.value,
        Assembly.GRCH38.value,
    )
    fusion_data_example_nonexonic = fusion_data_example_nonexonic(
        readData=ReadData(
            split=SplitReads(splitReads=100), spanning=SpanningReads(spanningReads=80)
        )
    )
    assert jaffa_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure
    assert jaffa_fusor_nonexonic.readData == fusion_data_example_nonexonic.readData

    # Test unknown partner
    jaffa.fusion_genes = "NA:PDGFRB"
    jaffa_fusor_unknown = await translator_instance.from_jaffa(
        jaffa, CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    assert jaffa_fusor_unknown.structure[0] == UnknownGeneElement()
    jaffa.fusion_genes = "TPM3:NA"
    jaffa_fusor_unknown = await translator_instance.from_jaffa(
        jaffa, CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    assert jaffa_fusor_unknown.structure[1] == UnknownGeneElement()


@pytest.mark.asyncio()
async def test_star_fusion(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test STAR-Fusion translator"""
    # Test exonic breakpoints
    star_fusion = STARFusion(
        left_gene="TPM3^ENSG00000143549.19",
        right_gene="PDGFRB^ENSG00000113721",
        left_breakpoint="chr1:154170465:-",
        right_breakpoint="chr5:150126612:-",
        annots='["INTERCHROMOSOMAL]',
        junction_read_count=100,
        spanning_frag_count=80,
    )

    star_fusion_fusor = await translator_instance.from_star_fusion(
        star_fusion,
        CoordinateType.INTER_RESIDUE.value,
        Assembly.GRCH38.value,
    )
    fusion_data_example = fusion_data_example(
        readData=ReadData(
            split=SplitReads(splitReads=100), spanning=SpanningReads(spanningReads=80)
        )
    )
    assert star_fusion_fusor.structure == fusion_data_example.structure
    assert star_fusion_fusor.readData == fusion_data_example.readData

    # Test non-exonic breakpoints
    star_fusion.left_breakpoint = "chr1:154173079:-"
    star_fusion.right_breakpoint = "chr5:150127173:-"

    star_fusion_fusor_nonexonic = await translator_instance.from_star_fusion(
        star_fusion,
        CoordinateType.RESIDUE.value,
        Assembly.GRCH38.value,
    )
    fusion_data_example_nonexonic = fusion_data_example_nonexonic(
        readData=ReadData(
            split=SplitReads(splitReads=100), spanning=SpanningReads(spanningReads=80)
        )
    )
    assert (
        star_fusion_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure
    )
    assert (
        star_fusion_fusor_nonexonic.readData == fusion_data_example_nonexonic.readData
    )

    # Test unknown partners
    star_fusion.left_gene = "NA"
    star_fusion_fusor_unknown = await translator_instance.from_star_fusion(
        star_fusion,
        CoordinateType.INTER_RESIDUE.value,
        Assembly.GRCH38.value,
    )
    assert star_fusion_fusor_unknown.structure[0] == UnknownGeneElement()
    star_fusion.left_gene = "TPM3"
    star_fusion.right_gene = "NA"
    star_fusion_fusor_unknown = await translator_instance.from_star_fusion(
        star_fusion,
        CoordinateType.INTER_RESIDUE.value,
        Assembly.GRCH38.value,
    )
    assert star_fusion_fusor_unknown.structure[1] == UnknownGeneElement()


@pytest.mark.asyncio()
async def test_fusion_catcher(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test Fusion Catcher translator"""
    # Test exonic breakpoint
    fusion_catcher = FusionCatcher(
        five_prime_partner="TPM3",
        three_prime_partner="PDGFRB",
        five_prime_fusion_point="1:154170465:-",
        three_prime_fusion_point="5:150126612:-",
        predicted_effect="exonic(no-known-CDS)/exonic(no-known-CDS)",
        spanning_unique_reads=100,
        spanning_reads=80,
        fusion_sequence="CTAGATGAC*TACTACTA",
    )

    fusion_catcher_fusor = await translator_instance.from_fusion_catcher(
        fusion_catcher,
        CoordinateType.INTER_RESIDUE.value,
        Assembly.GRCH38.value,
    )
    fusion_data_example = fusion_data_example(
        readData=ReadData(
            split=SplitReads(splitReads=100), spanning=SpanningReads(spanningReads=80)
        ),
        contig=ContigSequence(contig="CTAGATGAC*TACTACTA"),
    )
    assert fusion_catcher_fusor.structure == fusion_data_example.structure
    assert fusion_catcher_fusor.readData == fusion_data_example.readData
    assert fusion_catcher_fusor.contig == fusion_catcher_fusor.contig

    # Test non-exonic breakpoints
    fusion_catcher.five_prime_fusion_point = "1:154173079:-"
    fusion_catcher.three_prime_fusion_point = "5:150127173:-"

    fusion_catcher_fusor_nonexonic = await translator_instance.from_fusion_catcher(
        fusion_catcher,
        CoordinateType.RESIDUE.value,
        Assembly.GRCH38.value,
    )
    fusion_data_example_nonexonic = fusion_data_example_nonexonic(
        readData=ReadData(
            split=SplitReads(splitReads=100), spanning=SpanningReads(spanningReads=80)
        ),
        contig=ContigSequence(contig="CTAGATGAC*TACTACTA"),
    )
    assert (
        fusion_catcher_fusor_nonexonic.structure
        == fusion_data_example_nonexonic.structure
    )
    assert (
        fusion_catcher_fusor_nonexonic.readData
        == fusion_data_example_nonexonic.readData
    )
    assert fusion_catcher_fusor_nonexonic.contig == fusion_data_example_nonexonic.contig

    # Test unknown partners
    fusion_catcher.five_prime_partner = "NA"
    fusion_catcher_fusor_unknown = await translator_instance.from_fusion_catcher(
        fusion_catcher, CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    assert fusion_catcher_fusor_unknown.structure[0] == UnknownGeneElement()
    fusion_catcher.five_prime_partner = "TPM3"
    fusion_catcher.three_prime_partner = "NA"
    fusion_catcher_fusor_unknown = await translator_instance.from_fusion_catcher(
        fusion_catcher, CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    assert fusion_catcher_fusor_unknown.structure[1] == UnknownGeneElement()


@pytest.mark.asyncio()
async def test_fusion_map(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test Fusion Map translator"""
    # Test exonic breakpoint
    fusion_map_data = pl.DataFrame(
        {
            "KnownGene1": "TPM3",
            "KnownGene2": "PDGFRB",
            "Chromosome1": "1",
            "Position1": "154170465",
            "Chromosome2": "5",
            "Position2": "150126612",
            "FusionGene": "TPM3->PDGFRB",
            "SplicePatternClass": "CanonicalPattern[Major]",
            "FrameShiftClass": "InFrame",
        }
    )
    fusion_map_fusor = await translator_instance.from_fusion_map(
        fusion_map_data, CoordinateType.INTER_RESIDUE.value, Assembly.GRCH38.value
    )
    assert fusion_map_fusor.structure == fusion_data_example().structure

    # Test non-exonic breakpoint
    fusion_map_data_nonexonic = pl.DataFrame(
        {
            "KnownGene1": "TPM3",
            "KnownGene2": "PDGFRB",
            "Chromosome1": "1",
            "Position1": "154173079",
            "Chromosome2": "5",
            "Position2": "150127173",
            "FusionGene": "TPM3->PDGFRB",
            "SplicePatternClass": "CanonicalPattern[Major]",
            "FrameShiftClass": "InFrame",
        }
    )
    fusion_map_fusor_nonexonic = await translator_instance.from_fusion_map(
        fusion_map_data_nonexonic, CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    assert (
        fusion_map_fusor_nonexonic.structure
        == fusion_data_example_nonexonic().structure
    )


@pytest.mark.asyncio()
async def test_arriba(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test Arriba translator"""
    # Test exonic breakpoint
    arriba = Arriba(
        gene1="TPM3",
        gene2="PDGFRB",
        strand1="-/-",
        strand2="-/-",
        breakpoint1="1:154170465",
        breakpoint2="5:150126612",
        event_type="translocation",
        confidence="high",
        direction1="upstream",
        direction2="downstream",
        rf="in-frame",
        split_reads1=100,
        split_reads2=95,
        discordant_mates=30,
        coverage1=200,
        coverage2=190,
        fusion_transcript="CTAGATGAC_TACTACTA|GTACTACT",
    )

    arriba_fusor = await translator_instance.from_arriba(
        arriba,
        CoordinateType.INTER_RESIDUE.value,
        Assembly.GRCH38.value,
    )
    fusion_data_example = fusion_data_example(
        readData=ReadData(spanning=SpanningReads(spanningReads=30)),
        contig=ContigSequence(contig=arriba.fusion_transcript),
    )
    fusion_data_example.structure[0].coverage = BreakpointCoverage(fragmentCoverage=200)
    fusion_data_example.structure[0].anchoredReads = AnchoredReads(reads=100)
    fusion_data_example.structure[1].coverage = BreakpointCoverage(fragmentCoverage=190)
    fusion_data_example.structure[1].anchoredReads = AnchoredReads(reads=95)
    assert arriba_fusor.structure == fusion_data_example.structure
    assert arriba_fusor.readData == fusion_data_example.readData
    assert arriba_fusor.contig == fusion_data_example.contig

    # Test non-exonic breakpoint
    arriba.breakpoint1 = "1:154173079"
    arriba.breakpoint2 = "5:150127173"

    arriba_fusor_nonexonic = await translator_instance.from_arriba(
        arriba,
        CoordinateType.RESIDUE.value,
        Assembly.GRCH38.value,
    )
    fusion_data_example_nonexonic = fusion_data_example_nonexonic(
        readData=ReadData(spanning=SpanningReads(spanningReads=30)),
        contig=ContigSequence(contig=arriba.fusion_transcript),
    )
    fusion_data_example_nonexonic.structure[0].coverage = BreakpointCoverage(
        fragmentCoverage=200
    )
    fusion_data_example_nonexonic.structure[0].anchoredReads = AnchoredReads(reads=100)
    fusion_data_example_nonexonic.structure[1].coverage = BreakpointCoverage(
        fragmentCoverage=190
    )
    fusion_data_example_nonexonic.structure[1].anchoredReads = AnchoredReads(reads=95)
    assert arriba_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure
    assert arriba_fusor_nonexonic.readData == fusion_data_example_nonexonic.readData
    assert arriba_fusor_nonexonic.contig == fusion_data_example_nonexonic.contig

    # Test unknown partners
    arriba.gene1 = "NA"
    arriba_fusor_unknown = await translator_instance.from_arriba(
        arriba, CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    assert arriba_fusor_unknown.structure[0] == UnknownGeneElement()
    arriba.gene1 = "TPM3"
    arriba.gene2 = "NA"
    arriba_fusor_unknown = await translator_instance.from_arriba(
        arriba, CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    assert arriba_fusor_unknown.structure[1] == UnknownGeneElement()


@pytest.mark.asyncio()
async def test_cicero(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test CICERO translator"""
    # Test exonic breakpoint
    cicero = Cicero(
        gene_5prime="TPM3",
        gene_3prime="PDGFRB",
        chr_5prime="1",
        chr_3prime="5",
        pos_5prime=154170465,
        pos_3prime=150126612,
        sv_ort=">",
        event_type="CTX",
        reads_5prime=100,
        reads_3prime=90,
        coverage_5prime=200,
        coverage_3prime=190,
        contig="ATCATACTAGATACTACTACGATGAGAGAGTACATAGAT",
    )

    cicero_fusor = await translator_instance.from_cicero(
        cicero,
        CoordinateType.INTER_RESIDUE.value,
        Assembly.GRCH38.value,
    )
    fusion_data_example = fusion_data_example(
        contig=ContigSequence(contig=cicero.contig)
    )
    fusion_data_example.structure[0].coverage = BreakpointCoverage(fragmentCoverage=200)
    fusion_data_example.structure[0].anchoredReads = AnchoredReads(reads=100)
    fusion_data_example.structure[1].coverage = BreakpointCoverage(fragmentCoverage=190)
    fusion_data_example.structure[1].anchoredReads = AnchoredReads(reads=90)
    assert cicero_fusor.structure == fusion_data_example.structure
    assert cicero_fusor.readData == fusion_data_example.readData
    assert cicero_fusor.contig == fusion_data_example.contig

    # Test non-exonic breakpoint
    cicero.pos_5prime = 154173079
    cicero.pos_3prime = 150127173

    cicero_fusor_nonexonic = await translator_instance.from_cicero(
        cicero,
        CoordinateType.RESIDUE.value,
        Assembly.GRCH38.value,
    )
    fusion_data_example_nonexonic = fusion_data_example_nonexonic(
        contig=ContigSequence(contig=cicero.contig)
    )
    fusion_data_example_nonexonic.structure[0].coverage = BreakpointCoverage(
        fragmentCoverage=200
    )
    fusion_data_example_nonexonic.structure[0].anchoredReads = AnchoredReads(reads=100)
    fusion_data_example_nonexonic.structure[1].coverage = BreakpointCoverage(
        fragmentCoverage=190
    )
    fusion_data_example_nonexonic.structure[1].anchoredReads = AnchoredReads(reads=90)
    assert cicero_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure
    assert cicero_fusor_nonexonic.readData == fusion_data_example_nonexonic.readData
    assert cicero_fusor_nonexonic.contig == fusion_data_example_nonexonic.contig

    # Test case where the called fusion does not have confident biological meaning
    cicero.sv_ort = "?"

    non_confident_bio = await translator_instance.from_cicero(
        cicero,
        CoordinateType.RESIDUE.value,
        Assembly.GRCH38.value,
    )
    assert (
        non_confident_bio
        == "CICERO annotation indicates that this event does not have confident biological meaning"
    )

    # Test case where multiple gene symbols are reported for a fusion partner
    cicero.gene_3prime = "PDGFRB,PDGFRB-FGFR4,FGFR4"

    multiple_genes_fusion_partner = await translator_instance.from_cicero(
        cicero,
        CoordinateType.RESIDUE.value,
        Assembly.GRCH38.value,
    )
    assert (
        multiple_genes_fusion_partner
        == "Ambiguous gene symbols are reported by CICERO for at least one of the fusion partners"
    )

    # Test unknown partners
    cicero.sv_ort = ">"
    cicero.gene_5prime = "NA"
    cicero.gene_3prime = "PDGRB"
    cicero_fusor_unknown = await translator_instance.from_cicero(
        cicero, CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    assert cicero_fusor_unknown.structure[0] == UnknownGeneElement()
    cicero.gene_5prime = "TPM3"
    cicero.gene_3prime = "NA"
    cicero_fusor_unknown = await translator_instance.from_cicero(
        cicero, CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    assert cicero_fusor_unknown.structure[1] == UnknownGeneElement()


@pytest.mark.asyncio()
async def test_enfusion(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test EnFusion translator"""
    # Test exonic breakpoint
    enfusion = EnFusion(
        gene_5prime="TPM3",
        gene_3prime="PDGFRB",
        chr_5prime=1,
        chr_3prime=5,
        break_5prime=154170465,
        break_3prime=150126612,
    )

    enfusion_fusor = await translator_instance.from_enfusion(
        enfusion,
        CoordinateType.INTER_RESIDUE.value,
        Assembly.GRCH38.value,
    )
    assert enfusion_fusor.structure == fusion_data_example().structure

    # Test non-exonic breakpoint
    enfusion.break_5prime = 154173079
    enfusion.break_3prime = 150127173

    enfusion_fusor_nonexonic = await translator_instance.from_enfusion(
        enfusion,
        CoordinateType.RESIDUE.value,
        Assembly.GRCH38.value,
    )
    assert (
        enfusion_fusor_nonexonic.structure == fusion_data_example_nonexonic().structure
    )

    # Test unknown partner
    enfusion.gene_5prime = "NA"
    enfusion_fusor_unknown = await translator_instance.from_enfusion(
        enfusion, CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    assert enfusion_fusor_unknown.structure[0] == UnknownGeneElement()
    enfusion.gene_5prime = "TPM3"
    enfusion.gene_3prime = "NA"
    enfusion_fusor_unknown = await translator_instance.from_enfusion(
        enfusion, CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    assert enfusion_fusor_unknown.structure[1] == UnknownGeneElement()


@pytest.mark.asyncio()
async def test_genie(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test GENIE Translator"""
    # Test exonic breakpoint
    genie = Genie(
        site1_hugo="TPM3",
        site2_hugo="PDGFRB",
        site1_chrom=1,
        site2_chrom=5,
        site1_pos=154170465,
        site2_pos=150126612,
        annot="TMP3 (NM_152263.4) - PDGFRB (NM_002609.4) fusion",
        reading_frame="In_frame",
    )

    genie_fusor = await translator_instance.from_genie(
        genie,
        CoordinateType.INTER_RESIDUE.value,
        Assembly.GRCH38.value,
    )
    assert genie_fusor.structure == fusion_data_example().structure

    # Test non-exonic breakpoint
    genie.site1_pos = 154173079
    genie.site2_pos = 150127173

    genie_fusor_nonexonic = await translator_instance.from_genie(
        genie,
        CoordinateType.RESIDUE.value,
        Assembly.GRCH38.value,
    )
    assert genie_fusor_nonexonic.structure == fusion_data_example_nonexonic().structure

    # Test unknown partner
    genie.site1_hugo = "NA"
    genie_fusor_unknown = await translator_instance.from_genie(
        genie, CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    assert genie_fusor_unknown.structure[0] == UnknownGeneElement()
    genie.site1_hugo = "TPM3"
    genie.site2_hugo = "NA"
    genie_fusor_unknown = await translator_instance.from_genie(
        genie, CoordinateType.RESIDUE.value, Assembly.GRCH38.value
    )
    assert genie_fusor_unknown.structure[1] == UnknownGeneElement()
