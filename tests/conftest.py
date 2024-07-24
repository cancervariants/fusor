"""Module containing methods and fixtures used throughout tests."""

import logging

import pytest
from cool_seq_tool.app import CoolSeqTool

from fusor.fusor import FUSOR


def pytest_addoption(parser):
    """Add custom commands to pytest invocation.
    See https://docs.pytest.org/en/7.1.x/reference/reference.html#parser
    """
    parser.addoption(
        "--verbose-logs",
        action="store_true",
        default=False,
        help="show noisy module logs",
    )


def pytest_configure(config):
    """Configure pytest setup."""
    if not config.getoption("--verbose-logs"):
        logging.getLogger("botocore").setLevel(logging.INFO)
        logging.getLogger("boto3").setLevel(logging.INFO)
        logging.getLogger("urllib3").setLevel(logging.INFO)
        logging.getLogger("nose").setLevel(logging.INFO)


@pytest.fixture(scope="session")
def fusor_instance():
    """Create test fixture for fusor object

    Suppresses checks for CoolSeqTool external resources. Otherwise, on CST startup,
    it will try to check that its MANE summary file is up-to-date, which is an FTP call
    to the NCBI servers and can hang sometimes.

    If those files aren't available, create a CST instance in another session -- by
    default, it should save files to a centralized location that this test instance can
    access.
    """
    cst = CoolSeqTool(force_local_files=True)
    return FUSOR(cool_seq_tool=cst)


@pytest.fixture(scope="session")
def braf_gene():
    """Create gene params for BRAF."""
    return {
        "type": "Gene",
        "id": "hgnc:1097",
        "label": "BRAF",
        "mappings": [
            {
                "coding": {"code": "673", "system": "ncbigene"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "ENSG00000157764", "system": "ensembl"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "CCDS5863", "system": "ccds"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "1943", "system": "iuphar"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "119066", "system": "orphanet"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "BRAF", "system": "cosmic"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "2284096", "system": "pubmed"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "uc003vwc.5", "system": "ucsc"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "164757", "system": "omim"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "NM_004333", "system": "refseq"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "CCDS87555", "system": "ccds"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "P15056", "system": "uniprot"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "M95712", "system": "ena.embl"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "OTTHUMG00000157457", "system": "vega"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "1565476", "system": "pubmed"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "CCDS94219", "system": "ccds"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "CCDS94218", "system": "ccds"},
                "relation": "relatedMatch",
            },
        ],
        "alternativeLabels": ["BRAF1", "BRAF-1", "RAFB1", "NS7", "B-RAF1", "B-raf"],
        "extensions": [
            {
                "name": "approved_name",
                "value": "B-Raf proto-oncogene, serine/threonine kinase",
            },
            {
                "name": "ensembl_locations",
                "value": [
                    {
                        "type": "SequenceLocation",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        },
                        "start": 140719326,
                        "end": 140924929,
                    }
                ],
            },
            {
                "name": "ncbi_locations",
                "value": [
                    {
                        "type": "SequenceLocation",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        },
                        "start": 140713327,
                        "end": 140924929,
                    }
                ],
            },
            {"name": "ncbi_gene_type", "value": "protein-coding"},
            {
                "name": "hgnc_locus_type",
                "value": "gene with protein product",
            },
            {"name": "ensembl_biotype", "value": "protein_coding"},
            {"name": "strand", "value": "-"},
            {"name": "symbol_status", "value": "approved"},
        ],
    }


@pytest.fixture(scope="session")
def alk_gene():
    """Create test fixture for ALK gene params"""
    return {
        "id": "hgnc:427",
        "type": "Gene",
        "label": "ALK",
        "description": None,
        "alternativeLabels": ["NBLST3", "CD246", "ALK1"],
        "extensions": [
            {"name": "symbol_status", "value": "approved", "description": None},
            {
                "name": "approved_name",
                "value": "ALK receptor tyrosine kinase",
                "description": None,
            },
            {"name": "strand", "value": "-", "description": None},
            {
                "name": "ensembl_locations",
                "value": [
                    {
                        "id": "ga4gh:SL.V-yTsF-F4eHxeDHeU5KZIF3ZOzE2vUnG",
                        "type": "SequenceLocation",
                        "digest": "V-yTsF-F4eHxeDHeU5KZIF3ZOzE2vUnG",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
                        },
                        "start": 29192773,
                        "end": 29921586,
                    }
                ],
                "description": None,
            },
            {
                "name": "ncbi_locations",
                "value": [
                    {
                        "id": "ga4gh:SL.V-yTsF-F4eHxeDHeU5KZIF3ZOzE2vUnG",
                        "type": "SequenceLocation",
                        "digest": "V-yTsF-F4eHxeDHeU5KZIF3ZOzE2vUnG",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
                        },
                        "start": 29192773,
                        "end": 29921586,
                    }
                ],
                "description": None,
            },
            {
                "name": "hgnc_locus_type",
                "value": "gene with protein product",
                "description": None,
            },
            {"name": "ncbi_gene_type", "value": "protein-coding", "description": None},
            {"name": "ensembl_biotype", "value": "protein_coding", "description": None},
        ],
        "mappings": [
            {
                "coding": {
                    "label": None,
                    "system": "ensembl",
                    "version": None,
                    "code": "ENSG00000171094",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ncbigene",
                    "version": None,
                    "code": "238",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "orphanet",
                    "version": None,
                    "code": "160020",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "hcdmdb",
                    "version": None,
                    "code": "CD246",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ucsc",
                    "version": None,
                    "code": "uc002rmy.4",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "refseq",
                    "version": None,
                    "code": "NM_004304",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ccds",
                    "version": None,
                    "code": "CCDS33172",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "omim",
                    "version": None,
                    "code": "105590",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ena.embl",
                    "version": None,
                    "code": "D45915",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "vega",
                    "version": None,
                    "code": "OTTHUMG00000152034",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "uniprot",
                    "version": None,
                    "code": "Q9UM73",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "iuphar",
                    "version": None,
                    "code": "1839",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "pubmed",
                    "version": None,
                    "code": "8122112",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "cosmic",
                    "version": None,
                    "code": "ALK",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ccds",
                    "version": None,
                    "code": "CCDS86828",
                },
                "relation": "relatedMatch",
            },
        ],
    }


@pytest.fixture(scope="session")
def tpm3_gene():
    """Create test fixture for TPM3 gene"""
    return {
        "id": "hgnc:12012",
        "type": "Gene",
        "label": "TPM3",
        "description": None,
        "alternativeLabels": [
            "TM3",
            "NEM1~withdrawn",
            "TM30",
            "TM5",
            "TRK",
            "HEL-S-82p",
            "NEM1",
            "OK/SW-cl.5",
            "TM30nm",
            "hscp30",
            "FLJ35371",
            "TPMsk3",
            "HEL-189",
            "CFTD",
            "TPM3nu",
            "TM-5",
            "CAPM1",
        ],
        "extensions": [
            {"name": "symbol_status", "value": "approved", "description": None},
            {"name": "approved_name", "value": "tropomyosin 3", "description": None},
            {
                "name": "previous_symbols",
                "value": ["FLJ35371", "NEM1", "NEM1~withdrawn"],
                "description": None,
            },
            {"name": "strand", "value": "-", "description": None},
            {
                "name": "ensembl_locations",
                "value": [
                    {
                        "id": "ga4gh:SL.cgdnkG0tZq9SpwTHMWMG4sjT9JGXQ-Ap",
                        "type": "SequenceLocation",
                        "digest": "cgdnkG0tZq9SpwTHMWMG4sjT9JGXQ-Ap",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        },
                        "start": 154155307,
                        "end": 154194648,
                    }
                ],
                "description": None,
            },
            {
                "name": "ncbi_locations",
                "value": [
                    {
                        "id": "ga4gh:SL.aVsAgF9lwnjLgy-DXECiDgavt5F0OsYR",
                        "type": "SequenceLocation",
                        "digest": "aVsAgF9lwnjLgy-DXECiDgavt5F0OsYR",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        },
                        "start": 154155307,
                        "end": 154192100,
                    }
                ],
                "description": None,
            },
            {
                "name": "hgnc_locus_type",
                "value": "gene with protein product",
                "description": None,
            },
            {"name": "ncbi_gene_type", "value": "protein-coding", "description": None},
            {"name": "ensembl_biotype", "value": "protein_coding", "description": None},
        ],
        "mappings": [
            {
                "coding": {
                    "label": None,
                    "system": "ensembl",
                    "version": None,
                    "code": "ENSG00000143549",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ncbigene",
                    "version": None,
                    "code": "7170",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ccds",
                    "version": None,
                    "code": "CCDS41403",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ucsc",
                    "version": None,
                    "code": "uc001fec.3",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "pubmed",
                    "version": None,
                    "code": "25369766",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ccds",
                    "version": None,
                    "code": "CCDS41401",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ccds",
                    "version": None,
                    "code": "CCDS60275",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "uniprot",
                    "version": None,
                    "code": "P06753",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "cosmic",
                    "version": None,
                    "code": "TPM3",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ccds",
                    "version": None,
                    "code": "CCDS60274",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ena.embl",
                    "version": None,
                    "code": "BC008425",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ccds",
                    "version": None,
                    "code": "CCDS41402",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ccds",
                    "version": None,
                    "code": "CCDS1060",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "omim",
                    "version": None,
                    "code": "191030",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "orphanet",
                    "version": None,
                    "code": "120227",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ccds",
                    "version": None,
                    "code": "CCDS41400",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "refseq",
                    "version": None,
                    "code": "NM_152263",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "vega",
                    "version": None,
                    "code": "OTTHUMG00000035853",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "pubmed",
                    "version": None,
                    "code": "1829807",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ccds",
                    "version": None,
                    "code": "CCDS72922",
                },
                "relation": "relatedMatch",
            },
        ],
    }


@pytest.fixture(scope="module")
def exhaustive_example(alk_gene, braf_gene, tpm3_gene):
    """Create test fixture for a fake fusion exemplifying most major field types, in
    'expanded' form (ie properties augmented by VICC descriptors)
    """
    return {
        "type": "CategoricalFusion",
        "criticalFunctionalDomains": [
            {
                "type": "FunctionalDomain",
                "id": "interpro:IPR020635",
                "label": "Tyrosine-protein kinase, catalytic domain",
                "status": "lost",
                "associatedGene": alk_gene,
                "sequenceLocation": {
                    "id": "ga4gh:SL.aYx-iUOFEw7GVZb4fwrQLkQQahpiIAVp",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NP_004295.2",
                        "refgetAccession": "SQ.q9CnK-HKWh9eqhOi8FlzR7M0pCmUrWPs",
                        "type": "SequenceReference",
                    },
                    "start": 1116,
                    "end": 1383,
                },
            }
        ],
        "structure": [
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_152263.3",
                "exonStart": 1,
                "exonStartOffset": 0,
                "exonEnd": 8,
                "exonEndOffset": 0,
                "gene": tpm3_gene,
                "elementGenomicStart": {
                    "id": "ga4gh:SL.2K1vML0ofuYrYncrzzXUQOISRFJldZrO",
                    "description": None,
                    "xrefs": None,
                    "alternativeLabels": None,
                    "extensions": None,
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NC_000001.11",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        "type": "SequenceReference",
                    },
                    "start": 154192135,
                    "end": 154192136,
                },
                "elementGenomicEnd": {
                    "id": "ga4gh:SL.rtR6x2NnJEpROlxiT_DY9C-spf6ijYQi",
                    "description": None,
                    "xrefs": None,
                    "alternativeLabels": None,
                    "extensions": None,
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NC_000001.11",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        "type": "SequenceReference",
                    },
                    "start": 154170399,
                    "end": 154170400,
                },
            },
            {
                "type": "GeneElement",
                "gene": alk_gene,
            },
            {
                "type": "LinkerSequenceElement",
                "linkerSequence": {
                    "id": "fusor.sequence:ACGT",
                    "type": "LiteralSequenceExpression",
                    "label": None,
                    "description": None,
                    "xrefs": None,
                    "alternativeLabels": None,
                    "extensions": None,
                    "sequence_id": None,
                    "sequence": "ACGT",
                },
            },
            {
                "type": "TemplatedSequenceElement",
                "region": {
                    "id": "ga4gh:SL.gb3ew2XQ-Doi1AtvlmajeZO7fS1eDPg_",
                    "description": None,
                    "xrefs": None,
                    "alternativeLabels": None,
                    "extensions": None,
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NC_000023.11",
                        "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                        "type": "SequenceReference",
                    },
                    "start": 44908820,
                    "end": 44908822,
                },
                "strand": "+",
            },
            {"type": "MultiplePossibleGenesElement"},
        ],
        "regulatoryElement": {
            "type": "RegulatoryElement",
            "regulatoryClass": "promoter",
            "associatedGene": braf_gene,
        },
    }


@pytest.fixture()
def fusion_example():
    """Create test fixture for a fake fusion without additional property expansion."""
    return {
        "type": "CategoricalFusion",
        "readingFramePreserved": True,
        "criticalFunctionalDomains": [
            {
                "type": "FunctionalDomain",
                "id": "interpro:IPR020635",
                "label": "Tyrosine-protein kinase, catalytic domain",
                "status": "lost",
                "associatedGene": {
                    "id": "normalize.gene:hgnc%3A427",
                    "type": "Gene",
                    "label": "ALK",
                    "gene_id": "hgnc:427",
                },
                "sequenceLocation": {
                    "id": "ga4gh:SL.aYx-iUOFEw7GVZb4fwrQLkQQahpiIAVp",
                    "description": None,
                    "xrefs": None,
                    "alternativeLabels": None,
                    "extensions": None,
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NP_004295.2",
                        "refgetAccession": "SQ.q9CnK-HKWh9eqhOi8FlzR7M0pCmUrWPs",
                        "type": "SequenceReference",
                    },
                    "start": 1116,
                    "end": 1383,
                },
            }
        ],
        "structure": [
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_152263.3",
                "exonStart": 1,
                "exonStartOffset": 0,
                "exonEnd": 8,
                "exonEndOffset": 0,
                "gene": {
                    "type": "Gene",
                    "label": "TPM3",
                    "id": "hgnc:12012",
                },
                "elementGenomicStart": {
                    "id": "ga4gh:SL.2K1vML0ofuYrYncrzzXUQOISRFJldZrO",
                    "description": None,
                    "xrefs": None,
                    "alternativeLabels": None,
                    "extensions": None,
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NC_000001.11",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        "type": "SequenceReference",
                    },
                    "start": 154192135,
                    "end": 154192136,
                },
                "elementGenomicEnd": {
                    "id": "ga4gh:SL.rtR6x2NnJEpROlxiT_DY9C-spf6ijYQi",
                    "description": None,
                    "xrefs": None,
                    "alternativeLabels": None,
                    "extensions": None,
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NC_000001.11",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        "type": "SequenceReference",
                    },
                    "start": 154170399,
                    "end": 154170400,
                },
            },
            {
                "type": "GeneElement",
                "gene": {
                    "type": "Gene",
                    "label": "ALK",
                    "id": "hgnc:427",
                },
            },
            {
                "type": "LinkerSequenceElement",
                "linkerSequence": {
                    "id": "fusor.sequence:ACGT",
                    "type": "LiteralSequenceExpression",
                    "sequence": "ACGT",
                },
            },
            {
                "type": "TemplatedSequenceElement",
                "region": {
                    "id": "ga4gh:SL.gb3ew2XQ-Doi1AtvlmajeZO7fS1eDPg_",
                    "description": None,
                    "xrefs": None,
                    "alternativeLabels": None,
                    "extensions": None,
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NC_000023.11",
                        "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                        "type": "SequenceReference",
                    },
                    "start": 44908820,
                    "end": 44908822,
                },
                "strand": "+",
            },
            {"type": "MultiplePossibleGenesElement"},
        ],
        "regulatoryElement": {
            "type": "RegulatoryElement",
            "regulatoryClass": "promoter",
            "associatedGene": {
                "type": "Gene",
                "label": "BRAF",
                "id": "hgnc:1097",
            },
        },
    }
