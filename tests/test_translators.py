"""Module for testing FUSOR Translators"""

import polars as pl
import pytest
from cool_seq_tool.schemas import Assembly

from fusor.models import AssayedFusion
from fusor.translator import Caller


@pytest.fixture(scope="module")
def fusion_data_example():
    """Create example assayed fusion for TPM3::PDGFRB with exonic breakpoints"""
    params = {
        "type": "AssayedFusion",
        "structure": [
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_152263.4",
                "exonEnd": 8,
                "exonEndOffset": -66,
                "gene": {"id": "hgnc:12012", "type": "Gene", "label": "TPM3"},
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
                "gene": {"id": "hgnc:8804", "type": "Gene", "label": "PDGFRB"},
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
    return AssayedFusion(**params)


@pytest.fixture(scope="module")
def fusion_data_example_nonexonic():
    """Create example assayed fusion for TPM3::PDGFRB with non-exonic breakpoints"""
    params = {
        "type": "AssayedFusion",
        "structure": [
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_152263.4",
                "exonEnd": 4,
                "exonEndOffset": 5,
                "gene": {"id": "hgnc:12012", "type": "Gene", "label": "TPM3"},
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
                "gene": {"id": "hgnc:8804", "type": "Gene", "label": "PDGFRB"},
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
    return AssayedFusion(**params)


def test_gene_element_arriba(translator_instance):
    """Test gene selection for Arriba"""
    genes = "RP1-222H5.1(151985),MIR3672(13973)"
    gene = translator_instance._get_gene_element(genes=genes, caller=Caller.ARRIBA)
    assert gene.gene.label == "MIR3672"


@pytest.mark.asyncio()
async def test_jaffa(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test JAFFA translator"""
    # Test exonic breakpoint
    fusion_genes = "TPM3:PDGFRB"
    chrom1 = "chr1"
    base1 = 154170465
    chrom2 = "chr5"
    base2 = 150126612
    rearrangement = True
    classification = "HighConfidence"
    inframe = True

    jaffa_fusor = await translator_instance.from_jaffa(
        fusion_genes,
        chrom1,
        base1,
        chrom2,
        base2,
        rearrangement,
        classification,
        inframe,
        Assembly.GRCH38.value,
    )
    assert jaffa_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    fusion_genes = "TPM3:PDGFRB"
    chrom1 = "chr1"
    base1 = 154173078
    chrom2 = "chr5"
    base2 = 150127173
    rearrangement = True
    classification = "HighConfidence"
    inframe = True

    jaffa_fusor_nonexonic = await translator_instance.from_jaffa(
        fusion_genes,
        chrom1,
        base1,
        chrom2,
        base2,
        rearrangement,
        classification,
        inframe,
        Assembly.GRCH38.value,
    )
    assert jaffa_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure


@pytest.mark.asyncio()
async def test_star_fusion(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test STAR-Fusion translator"""
    # Test exonic breakpoints
    left_gene = "TPM3^ENSG00000143549.19"
    right_gene = "PDGFRB^ENSG00000113721"
    left_breakpoint = "chr1:154170465:-"
    right_breakpoint = "chr5:150126612:-"
    annots = '["INTERCHROMOSOMAL]'

    star_fusion_fusor = await translator_instance.from_star_fusion(
        left_gene,
        right_gene,
        left_breakpoint,
        right_breakpoint,
        annots,
        Assembly.GRCH38.value,
    )
    assert star_fusion_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoints
    left_gene = "TPM3^ENSG00000143549.19"
    right_gene = "PDGFRB^ENSG00000113721"
    left_breakpoint = "chr1:154173078:-"
    right_breakpoint = "chr5:150127173:-"
    annots = '["INTERCHROMOSOMAL]'

    star_fusion_fusor_nonexonic = await translator_instance.from_star_fusion(
        left_gene,
        right_gene,
        left_breakpoint,
        right_breakpoint,
        annots,
        Assembly.GRCH38.value,
    )
    assert (
        star_fusion_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure
    )


@pytest.mark.asyncio()
async def test_fusion_catcher(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test Fusion Catcher translator"""
    # Test exonic breakpoint
    five_prime_partner = "TPM3"
    three_prime_partner = "PDGFRB"
    five_prime_fusion_point = "1:154170465:-"
    three_prime_fusion_point = "5:150126612:-"
    predicted_effect = "exonic(no-known-CDS)/exonic(no-known-CDS)"

    fusion_catcher_fusor = await translator_instance.from_fusion_catcher(
        five_prime_partner,
        three_prime_partner,
        five_prime_fusion_point,
        three_prime_fusion_point,
        predicted_effect,
        Assembly.GRCH38.value,
    )
    assert fusion_catcher_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    five_prime_partner = "TPM3"
    three_prime_partner = "PDGFRB"
    five_prime_fusion_point = "1:154173078:-"
    three_prime_fusion_point = "5:150127173:-"
    predicted_effect = "exonic(no-known-CDS)/exonic(no-known-CDS)"

    fusion_catcher_fusor_nonexonic = await translator_instance.from_fusion_catcher(
        five_prime_partner,
        three_prime_partner,
        five_prime_fusion_point,
        three_prime_fusion_point,
        predicted_effect,
        Assembly.GRCH38.value,
    )
    assert (
        fusion_catcher_fusor_nonexonic.structure
        == fusion_data_example_nonexonic.structure
    )


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
        fusion_map_data, Assembly.GRCH38.value
    )
    assert fusion_map_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    fusion_map_data_nonexonic = pl.DataFrame(
        {
            "KnownGene1": "TPM3",
            "KnownGene2": "PDGFRB",
            "Chromosome1": "1",
            "Position1": "154173078",
            "Chromosome2": "5",
            "Position2": "150127173",
            "FusionGene": "TPM3->PDGFRB",
            "SplicePatternClass": "CanonicalPattern[Major]",
            "FrameShiftClass": "InFrame",
        }
    )
    fusion_map_fusor_nonexonic = await translator_instance.from_fusion_map(
        fusion_map_data_nonexonic, Assembly.GRCH38.value
    )
    assert (
        fusion_map_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure
    )


@pytest.mark.asyncio()
async def test_arriba(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test Arriba translator"""
    # Test exonic breakpoint
    gene1 = "TPM3"
    gene2 = "PDGFRB"
    strand1 = "-/-"
    strand2 = "-/-"
    breakpoint1 = "1:154170465"
    breakpoint2 = "5:150126612"
    event = "translocation"
    confidence = "high"
    direction1 = "upstream"
    direction2 = "downstream"
    rf = "in-frame"

    arriba_fusor = await translator_instance.from_arriba(
        gene1,
        gene2,
        strand1,
        strand2,
        breakpoint1,
        breakpoint2,
        event,
        confidence,
        direction1,
        direction2,
        rf,
        Assembly.GRCH38.value,
    )
    assert arriba_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    gene1 = "TPM3"
    gene2 = "PDGFRB"
    strand1 = "-/-"
    strand2 = "-/-"
    breakpoint1 = "1:154173078"
    breakpoint2 = "5:150127173"
    event = "translocation"
    confidence = "high"
    direction1 = "upstream"
    direction2 = "downstream"
    rf = "in-frame"

    arriba_fusor_nonexonic = await translator_instance.from_arriba(
        gene1,
        gene2,
        strand1,
        strand2,
        breakpoint1,
        breakpoint2,
        event,
        confidence,
        direction1,
        direction2,
        rf,
        Assembly.GRCH38.value,
    )
    assert arriba_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure


@pytest.mark.asyncio()
async def test_cicero(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test CICERO translator"""
    # Test exonic breakpoint
    gene_5prime = "TPM3"
    gene_3prime = "PDGFRB"
    chr_5prime = "1"
    chr_3prime = "5"
    pos_5prime = 154170465
    pos_3prime = 150126612
    sv_ort = ">"
    event_type = "CTX"

    cicero_fusor = await translator_instance.from_cicero(
        gene_5prime,
        gene_3prime,
        chr_5prime,
        chr_3prime,
        pos_5prime,
        pos_3prime,
        sv_ort,
        event_type,
        Assembly.GRCH38.value,
    )
    assert cicero_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    gene_5prime = "TPM3"
    gene_3prime = "PDGFRB"
    chr_5prime = "1"
    chr_3prime = "5"
    pos_5prime = 154173078
    pos_3prime = 150127173
    sv_ort = ">"
    event_type = "CTX"

    cicero_fusor_nonexonic = await translator_instance.from_cicero(
        gene_5prime,
        gene_3prime,
        chr_5prime,
        chr_3prime,
        pos_5prime,
        pos_3prime,
        sv_ort,
        event_type,
        Assembly.GRCH38.value,
    )
    assert cicero_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure

    # Test case where the called fusion does not have confident biological meaning
    gene_5prime = "TPM3"
    gene_3prime = "PDGFRB"
    chr_5prime = "1"
    chr_3prime = "5"
    pos_5prime = 154173078
    pos_3prime = 150127173
    sv_ort = "?"
    event_type = "CTX"

    non_confident_bio = await translator_instance.from_cicero(
        gene_5prime,
        gene_3prime,
        chr_5prime,
        chr_3prime,
        pos_5prime,
        pos_3prime,
        sv_ort,
        event_type,
        Assembly.GRCH38.value,
    )
    assert (
        non_confident_bio
        == "CICERO annotation indicates that this event does not have confident biological meaning"
    )

    # Test case where multiple gene symbols are reported for a fusion partner
    gene_5prime = "TPM3"
    gene_3prime = "PDGFRB,PDGFRB-FGFR4,FGFR4"
    chr_5prime = "1"
    chr_3prime = "5"
    pos_5prime = 154173078
    pos_3prime = 150127173
    sv_ort = "?"
    event_type = "CTX"

    multiple_genes_fusion_partner = await translator_instance.from_cicero(
        gene_5prime,
        gene_3prime,
        chr_5prime,
        chr_3prime,
        pos_5prime,
        pos_3prime,
        sv_ort,
        event_type,
        Assembly.GRCH38.value,
    )
    assert (
        multiple_genes_fusion_partner
        == "Ambiguous gene symbols are reported by CICERO for at least one of the fusion partners"
    )


@pytest.mark.asyncio()
async def test_enfusion(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test EnFusion translator"""
    # Test exonic breakpoint
    gene_5prime = "TPM3"
    gene_3prime = "PDGFRB"
    chr_5prime = 1
    chr_3prime = 5
    break_5prime = 154170465
    break_3prime = 150126612

    enfusion_fusor = await translator_instance.from_enfusion(
        gene_5prime,
        gene_3prime,
        chr_5prime,
        chr_3prime,
        break_5prime,
        break_3prime,
        Assembly.GRCH38.value,
    )
    assert enfusion_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    gene_5prime = "TPM3"
    gene_3prime = "PDGFRB"
    chr_5prime = 1
    chr_3prime = 5
    break_5prime = 154173078
    break_3prime = 150127173

    enfusion_fusor_nonexonic = await translator_instance.from_enfusion(
        gene_5prime,
        gene_3prime,
        chr_5prime,
        chr_3prime,
        break_5prime,
        break_3prime,
        Assembly.GRCH38.value,
    )
    assert enfusion_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure


@pytest.mark.asyncio()
async def test_genie(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test GENIE Translator"""
    # Test exonic breakpoint
    site1_hugo = "TPM3"
    site2_hugo = "PDGFRB"
    site1_chrom = 1
    site2_chrom = 5
    site1_pos = 154170465
    site2_pos = 150126612
    annot = "TMP3 (NM_152263.4) - PDGFRB (NM_002609.4) fusion"
    reading_frame = "In_frame"

    genie_fusor = await translator_instance.from_genie(
        site1_hugo,
        site2_hugo,
        site1_chrom,
        site2_chrom,
        site1_pos,
        site2_pos,
        annot,
        reading_frame,
        Assembly.GRCH38.value,
    )
    assert genie_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    site1_hugo = "TPM3"
    site2_hugo = "PDGFRB"
    site1_chrom = 1
    site2_chrom = 5
    site1_pos = 154173078
    site2_pos = 150127173
    annot = "TMP3 (NM_152263.4) - PDGFRB (NM_002609.4) fusion"
    reading_frame = "In_frame"

    genie_fusor_nonexonic = await translator_instance.from_genie(
        site1_hugo,
        site2_hugo,
        site1_chrom,
        site2_chrom,
        site1_pos,
        site2_pos,
        annot,
        reading_frame,
        Assembly.GRCH38.value,
    )
    assert genie_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure
