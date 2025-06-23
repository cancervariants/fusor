"""Module containing methods and fixtures used throughout tests."""

import logging
from pathlib import Path

import pytest
from cool_seq_tool.app import CoolSeqTool

from fusor.fusor import FUSOR
from fusor.translator import Translator

FIXTURE_DATA_DIR = Path(__file__).parents[0].resolve() / "fixtures"


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
def fixture_data_dir():
    """Provide test data directory."""
    return FIXTURE_DATA_DIR


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


@pytest.fixture()
def translator_instance():
    """Create test fixture for translator object"""
    return Translator(fusor=FUSOR())


@pytest.fixture(scope="session")
def braf_gene():
    """Create gene params for BRAF."""
    return {
        "primaryCoding": {
            "id": "hgnc:1097",
            "code": "HGNC:1097",
            "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
        },
        "name": "BRAF",
        "conceptType": "Gene",
        "extensions": [
            {
                "name": "aliases",
                "value": ["B-RAF1", "NS7", "BRAF-1", "RAFB1", "BRAF1", "B-raf"],
            },
            {"name": "symbol_status", "value": "approved", "description": None},
            {
                "name": "approved_name",
                "value": "B-Raf proto-oncogene, serine/threonine kinase",
                "description": None,
            },
            {"name": "strand", "value": "-", "description": None},
            {
                "name": "ensembl_locations",
                "value": [
                    {
                        "id": "ga4gh:SL.fUv91vYrVHBMg-B_QW7UpOQj50g_49hb",
                        "type": "SequenceLocation",
                        "digest": "fUv91vYrVHBMg-B_QW7UpOQj50g_49hb",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        },
                        "start": 140719326,
                        "end": 140924929,
                    }
                ],
                "description": None,
            },
            {
                "name": "ncbi_locations",
                "value": [
                    {
                        "id": "ga4gh:SL.0nPwKHYNnTmJ06G-gSmz8BEhB_NTp-0B",
                        "type": "SequenceLocation",
                        "digest": "0nPwKHYNnTmJ06G-gSmz8BEhB_NTp-0B",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        },
                        "start": 140713327,
                        "end": 140924929,
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
                    "system": "ensembl",
                    "code": "ENSG00000157764",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ncbigene",
                    "code": "673",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "cosmic",
                    "code": "BRAF",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ena.embl",
                    "code": "M95712",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "omim",
                    "code": "164757",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "iuphar",
                    "code": "1943",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ucsc",
                    "code": "uc003vwc.5",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "vega",
                    "code": "OTTHUMG00000157457",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ccds",
                    "code": "CCDS87555",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "uniprot",
                    "code": "P15056",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "refseq",
                    "code": "NM_004333",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "pubmed",
                    "code": "1565476",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "orphanet",
                    "code": "119066",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "pubmed",
                    "code": "2284096",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ccds",
                    "code": "CCDS5863",
                },
                "relation": "relatedMatch",
            },
        ],
    }


@pytest.fixture(scope="session")
def alk_gene():
    """Create test fixture for ALK gene params"""
    return {
        "primaryCoding": {
            "id": "hgnc:427",
            "code": "HGNC:427",
            "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
        },
        "conceptType": "Gene",
        "name": "ALK",
        "extensions": [
            {
                "name": "aliases",
                "value": [
                    "WAP4",
                    "MPI",
                    "ALK1",
                    "ALP",
                    "HUSI-I",
                    "BLPI",
                    "WFDC4",
                    "HUSI",
                ],
            },
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
                    "system": "ensembl",
                    "code": "ENSG00000171094",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ncbigene",
                    "code": "238",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "orphanet",
                    "code": "160020",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "hcdmdb",
                    "code": "CD246",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ucsc",
                    "code": "uc002rmy.4",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "refseq",
                    "code": "NM_004304",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ccds",
                    "code": "CCDS33172",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "omim",
                    "code": "105590",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ena.embl",
                    "code": "D45915",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "vega",
                    "code": "OTTHUMG00000152034",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "uniprot",
                    "code": "Q9UM73",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "iuphar",
                    "code": "1839",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "pubmed",
                    "code": "8122112",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "cosmic",
                    "code": "ALK",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ccds",
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
        "primaryCoding": {
            "id": "hgnc:12012",
            "code": "HGNC:12012",
            "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
        },
        "conceptType": "Gene",
        "name": "TPM3",
        "extensions": [
            {
                "name": "aliases",
                "value": [
                    "TM30",
                    "CAPM1",
                    "HEL-S-82p",
                    "NEM1",
                    "TRK",
                    "CMYP4B",
                    "HEL-189",
                    "TM5",
                    "CFTD",
                    "TPMsk3",
                    "CMYP4A",
                    "TM3",
                    "hscp30",
                    "OK/SW-cl.5",
                    "NEM1~withdrawn",
                    "TM-5",
                    "CMYO4B",
                    "FLJ35371",
                    "TPM3nu",
                    "TM30nm",
                    "CMYO4A",
                ],
            },
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
                    "system": "ensembl",
                    "code": "ENSG00000143549",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ncbigene",
                    "code": "7170",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ccds",
                    "code": "CCDS41403",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ucsc",
                    "code": "uc001fec.3",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "pubmed",
                    "code": "25369766",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ccds",
                    "code": "CCDS41401",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ccds",
                    "code": "CCDS60275",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "uniprot",
                    "code": "P06753",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "cosmic",
                    "code": "TPM3",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ccds",
                    "code": "CCDS60274",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ena.embl",
                    "code": "BC008425",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ccds",
                    "code": "CCDS41402",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ccds",
                    "code": "CCDS1060",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "omim",
                    "code": "191030",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "orphanet",
                    "code": "120227",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ccds",
                    "code": "CCDS41400",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "refseq",
                    "code": "NM_152263",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "vega",
                    "code": "OTTHUMG00000035853",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "pubmed",
                    "code": "1829807",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "system": "ccds",
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
                "strand": -1,
                "exonStart": 1,
                "exonStartOffset": 0,
                "exonEnd": 8,
                "exonEndOffset": 0,
                "gene": tpm3_gene,
                "elementGenomicStart": {
                    "id": "ga4gh:SL.Q8vkGp7_xR9vI0PQ7g1IvUUeQ4JlJG8l",
                    "digest": "Q8vkGp7_xR9vI0PQ7g1IvUUeQ4JlJG8l",
                    "description": None,
                    "extensions": None,
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NC_000001.11",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        "type": "SequenceReference",
                    },
                    "end": 154192135,
                },
                "elementGenomicEnd": {
                    "id": "ga4gh:SL.Lnne0bSsgjzmNkKsNnXg98FeJSrDJuLb",
                    "digest": "Lnne0bSsgjzmNkKsNnXg98FeJSrDJuLb",
                    "description": None,
                    "extensions": None,
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NC_000001.11",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        "type": "SequenceReference",
                    },
                    "start": 154170399,
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
                    "description": None,
                    "extensions": None,
                    "sequence": "ACGT",
                },
            },
            {
                "type": "TemplatedSequenceElement",
                "region": {
                    "id": "ga4gh:SL.gb3ew2XQ-Doi1AtvlmajeZO7fS1eDPg_",
                    "description": None,
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
                "strand": 1,
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
                    "conceptType": "Gene",
                    "name": "ALK",
                    "primaryCoding": {
                        "id": "hgnc:427",
                        "code": "HGNC:427",
                        "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                    },
                },
                "sequenceLocation": {
                    "id": "ga4gh:SL.aYx-iUOFEw7GVZb4fwrQLkQQahpiIAVp",
                    "description": None,
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
                "strand": -1,
                "exonStart": 1,
                "exonStartOffset": 0,
                "exonEnd": 8,
                "exonEndOffset": 0,
                "gene": {
                    "conceptType": "Gene",
                    "name": "TPM3",
                    "primaryCoding": {
                        "id": "hgnc:12012",
                        "code": "HGNC:12012",
                        "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                    },
                },
                "elementGenomicStart": {
                    "id": "ga4gh:SL.Q8vkGp7_xR9vI0PQ7g1IvUUeQ4JlJG8l",
                    "digest": "Q8vkGp7_xR9vI0PQ7g1IvUUeQ4JlJG8l",
                    "description": None,
                    "extensions": None,
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NC_000001.11",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        "type": "SequenceReference",
                    },
                    "end": 154192135,
                },
                "elementGenomicEnd": {
                    "id": "ga4gh:SL.Lnne0bSsgjzmNkKsNnXg98FeJSrDJuLb",
                    "digest": "Lnne0bSsgjzmNkKsNnXg98FeJSrDJuLb",
                    "description": None,
                    "extensions": None,
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NC_000001.11",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        "type": "SequenceReference",
                    },
                    "start": 154170399,
                },
            },
            {
                "type": "GeneElement",
                "gene": {
                    "conceptType": "Gene",
                    "name": "ALK",
                    "primaryCoding": {
                        "id": "hgnc:427",
                        "code": "HGNC:427",
                        "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                    },
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
                "strand": 1,
            },
            {"type": "MultiplePossibleGenesElement"},
        ],
        "regulatoryElement": {
            "type": "RegulatoryElement",
            "regulatoryClass": "promoter",
            "associatedGene": {
                "conceptType": "Gene",
                "name": "BRAF",
                "primaryCoding": {
                    "id": "hgnc:1097",
                    "code": "HGNC:1097",
                    "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                },
            },
        },
    }
