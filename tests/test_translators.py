"""Module for testing FUSOR Translators"""

import polars as pl
import pytest
from cool_seq_tool.schemas import Assembly

from fusor.models import AssayedFusion


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
    gene = translator_instance._get_gene_element(genelist=genes, caller="arriba")
    assert gene[0].gene.label == "MIR3672"


@pytest.mark.asyncio()
async def test_jaffa(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test JAFFA translator"""
    # Test exonic breakpoint
    jaffa_data = pl.DataFrame(
        {
            "fusion genes": "TPM3:PDGFRB",
            "chrom1": "chr1",
            "base1": "154170465",
            "chrom2": "chr5",
            "base2": "150126612",
            "rearrangement": "TRUE",
            "classification": "HighConfidence",
            "inframe": "TRUE",
        }
    )

    jaffa_fusor = (
        await translator_instance.from_jaffa(jaffa_data, Assembly("GRCh38"))
    )[0]
    assert jaffa_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    jaffa_data_nonexonic = pl.DataFrame(
        {
            "fusion genes": "TPM3:PDGFRB",
            "chrom1": "chr1",
            "base1": "154173078",
            "chrom2": "chr5",
            "base2": "150127173",
            "rearrangement": "TRUE",
            "classification": "HighConfidence",
            "inframe": "TRUE",
        }
    )
    jaffa_fusor_nonexonic = (
        await translator_instance.from_jaffa(jaffa_data_nonexonic, Assembly("GRCh38"))
    )[0]
    assert jaffa_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure


@pytest.mark.asyncio()
async def test_star_fusion(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test STAR-Fusion translator"""
    # Test exonic breakpoints
    star_fusion_data = pl.DataFrame(
        {
            "LeftGene": "TPM3^ENSG00000143549.19",
            "RightGene": "PDGFRB^ENSG00000113721",
            "LeftBreakpoint": "chr1:154170465:-",
            "RightBreakpoint": "chr5:150126612:-",
            "annots": '["INTRACHROMOSOMAL[chr16:0.23Mb]"]',
        }
    )
    star_fusion_fusor = (
        await translator_instance.from_star_fusion(star_fusion_data, Assembly("GRCh38"))
    )[0]
    assert star_fusion_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoints
    star_fusion_data_nonexonic = pl.DataFrame(
        {
            "LeftGene": "TPM3^ENSG00000143549.19",
            "RightGene": "PDGFRB^ENSG00000113721",
            "LeftBreakpoint": "chr1:154173078:-",
            "RightBreakpoint": "chr5:150127173:-",
            "annots": '["INTRACHROMOSOMAL[chr16:0.23Mb]"]',
        }
    )
    star_fusion_fusor_nonexonic = (
        await translator_instance.from_star_fusion(
            star_fusion_data_nonexonic, Assembly("GRCh38")
        )
    )[0]
    assert (
        star_fusion_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure
    )


@pytest.mark.asyncio()
async def test_fusion_catcher(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test Fusion Catcher translator"""
    # Test exonic breakpoint
    fusion_catcher_data = pl.DataFrame(
        {
            "Gene_1_symbol(5end_fusion_partner)": "TPM3",
            "Gene_2_symbol(3end_fusion_partner)": "PDGFRB",
            "Fusion_point_for_gene_1(5end_fusion_partner)": "1:154170465:-",
            "Fusion_point_for_gene_2(3end_fusion_partner)": "5:150126612:-",
            "Predicted_effect": "exonic(no-known-CDS)/exonic(no-known-CDS",
        }
    )
    fusion_catcher_fusor = (
        await translator_instance.from_fusion_catcher(
            fusion_catcher_data, Assembly("GRCh38")
        )
    )[0]
    assert fusion_catcher_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    fusion_catcher_data_nonexonic = pl.DataFrame(
        {
            "Gene_1_symbol(5end_fusion_partner)": "TPM3",
            "Gene_2_symbol(3end_fusion_partner)": "PDGFRB",
            "Fusion_point_for_gene_1(5end_fusion_partner)": "1:154173078:-",
            "Fusion_point_for_gene_2(3end_fusion_partner)": "5:150127173:-",
            "Predicted_effect": "exonic(no-known-CDS)/exonic(no-known-CDS",
        }
    )
    fusion_catcher_fusor_nonexonic = (
        await translator_instance.from_fusion_catcher(
            fusion_catcher_data_nonexonic, Assembly("GRCh38")
        )
    )[0]
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
    fusion_map_fusor = (
        await translator_instance.from_fusion_map(fusion_map_data, Assembly("GRCh38"))
    )[0]
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
    fusion_map_fusor_nonexonic = (
        await translator_instance.from_fusion_map(
            fusion_map_data_nonexonic, Assembly("GRCh38")
        )
    )[0]
    assert (
        fusion_map_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure
    )


@pytest.mark.asyncio()
async def test_arriba(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test Arriba translator"""
    # Test exonic breakpoint
    arriba_data = pl.DataFrame(
        {
            "#gene1": "TPM3",
            "gene2": "PDGFRB",
            "breakpoint1": "1:154170465",
            "breakpoint2": "5:150126612",
            "type": ".",
            "confidence": "high",
            "reading_frame": "in-frame",
        }
    )
    arriba_fusor = (
        await translator_instance.from_arriba(arriba_data, Assembly("GRCh38"))
    )[0]
    assert arriba_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    arriba_data_nonexonic = pl.DataFrame(
        {
            "#gene1": "TPM3",
            "gene2": "PDGFRB",
            "breakpoint1": "1:154173078",
            "breakpoint2": "5:150127173",
            "type": ".",
            "confidence": "high",
            "reading_frame": "in-frame",
        }
    )
    arriba_fusor_nonexonic = (
        await translator_instance.from_arriba(arriba_data_nonexonic, Assembly("GRCh38"))
    )[0]
    assert arriba_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure


@pytest.mark.asyncio()
async def test_cicero(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test CICERO translator"""
    # Test exonic breakpoint
    cicero_data = pl.DataFrame(
        {
            "geneA": "TPM3",
            "geneB": "PDGFRB",
            "chrA": "1",
            "posA": "154170465",
            "chrB": "5",
            "posB": "150126612",
            "type": "CTX",
        }
    )
    cicero_fusor = (
        await translator_instance.from_cicero(cicero_data, Assembly("GRCh38"))
    )[0]
    assert cicero_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    cicero_data_nonexonic = pl.DataFrame(
        {
            "geneA": "TPM3",
            "geneB": "PDGFRB",
            "chrA": "1",
            "posA": "154173078",
            "chrB": "5",
            "posB": "150127173",
            "type": "CTX",
        }
    )
    cicero_fusor_nonexonic = (
        await translator_instance.from_cicero(cicero_data_nonexonic, Assembly("GRCh38"))
    )[0]
    assert cicero_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure


@pytest.mark.asyncio()
async def test_enfusion(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test EnFusion translator"""
    # Test exonic breakpoint
    enfusion_data = pl.DataFrame(
        {
            "Gene1": "TPM3",
            "Gene2": "PDGFRB",
            "Chr1": "1",
            "Break1": "154170465",
            "Chr2": "5",
            "Break2": "150126612",
        }
    )
    enfusion_fusor = (
        await translator_instance.from_enfusion(enfusion_data, Assembly("GRCh38"))
    )[0]
    assert enfusion_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    enfusion_data_nonexonic = pl.DataFrame(
        {
            "Gene1": "TPM3",
            "Gene2": "PDGFRB",
            "Chr1": "1",
            "Break1": "154173078",
            "Chr2": "5",
            "Break2": "150127173",
        }
    )
    enfusion_fusor_nonexonic = (
        await translator_instance.from_enfusion(
            enfusion_data_nonexonic, Assembly("GRCh38")
        )
    )[0]
    assert enfusion_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure


@pytest.mark.asyncio()
async def test_genie(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test GENIE Translator"""
    # Test exonic breakpoint
    genie_data = pl.DataFrame(
        {
            "Site1_Hugo_Symbol": "TPM3",
            "Site2_Hugo_Symbol": "PDGFRB",
            "Site1_Chromosome": "1",
            "Site1_Position": "154170465",
            "Site2_Chromosome": "5",
            "Site2_Position": "150126612",
            "Annotation": "TMP3 (NM_152263.4) - PDGFRB (NM_002609.4) fusion",
            "Site2_Effect_on_Frame": "In_frame",
        }
    )
    genie_fusor = (
        await translator_instance.from_genie(genie_data, Assembly("GRCh38"))
    )[0]
    assert genie_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    genie_data_nonexonic = pl.DataFrame(
        {
            "Site1_Hugo_Symbol": "TPM3",
            "Site2_Hugo_Symbol": "PDGFRB",
            "Site1_Chromosome": "1",
            "Site1_Position": "154173078",
            "Site2_Chromosome": "5",
            "Site2_Position": "150127173",
            "Annotation": "TMP3 (NM_152263.4) - PDGFRB (NM_002609.4) fusion",
            "Site2_Effect_on_Frame": "In_frame",
        }
    )
    genie_fusor_nonexonic = (
        await translator_instance.from_genie(genie_data_nonexonic, Assembly("GRCh38"))
    )[0]
    assert genie_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure
