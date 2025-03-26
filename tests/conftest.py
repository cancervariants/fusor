"""Module containing methods and fixtures used throughout tests."""

import asyncio
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
def event_loop():
    """Create an instance of the default event loop for each test case."""
    loop = asyncio.get_event_loop_policy().new_event_loop()
    yield loop
    loop.close()


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


@pytest.fixture(scope="session")
def translator_instance():
    """Create test fixture for translator object"""
    return Translator(fusor=FUSOR())


@pytest.fixture(scope="session")
def braf_gene():
    """Create gene params for BRAF."""
    return {
        "concept_id": "hgnc:1097",
        "type": "Gene",
        "symbol": "BRAF",
        "description": None,
        "alternativeLabels": ["NS7", "BRAF1", "RAFB1", "B-RAF1", "BRAF-1", "B-raf"],
        "extensions": [
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
                    "label": None,
                    "system": "ensembl",
                    "version": None,
                    "code": "ENSG00000157764",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ncbigene",
                    "version": None,
                    "code": "673",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "cosmic",
                    "version": None,
                    "code": "BRAF",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ena.embl",
                    "version": None,
                    "code": "M95712",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "omim",
                    "version": None,
                    "code": "164757",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "iuphar",
                    "version": None,
                    "code": "1943",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ucsc",
                    "version": None,
                    "code": "uc003vwc.5",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "vega",
                    "version": None,
                    "code": "OTTHUMG00000157457",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ccds",
                    "version": None,
                    "code": "CCDS87555",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "uniprot",
                    "version": None,
                    "code": "P15056",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "refseq",
                    "version": None,
                    "code": "NM_004333",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "pubmed",
                    "version": None,
                    "code": "1565476",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "orphanet",
                    "version": None,
                    "code": "119066",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "pubmed",
                    "version": None,
                    "code": "2284096",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ccds",
                    "version": None,
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
        "concept_id": "hgnc:427",
        "type": "Gene",
        "symbol": "ALK",
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
        "concept_id": "hgnc:12012",
        "type": "Gene",
        "symbol": "TPM3",
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
                    "id": "ga4gh:SL.Q8vkGp7_xR9vI0PQ7g1IvUUeQ4JlJG8l",
                    "digest": "Q8vkGp7_xR9vI0PQ7g1IvUUeQ4JlJG8l",
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
                    "end": 154192135,
                },
                "elementGenomicEnd": {
                    "id": "ga4gh:SL.Lnne0bSsgjzmNkKsNnXg98FeJSrDJuLb",
                    "digest": "Lnne0bSsgjzmNkKsNnXg98FeJSrDJuLb",
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
                    "type": "Gene",
                    "symbol": "ALK",
                    "concept_id": "hgnc:427",
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
                    "symbol": "TPM3",
                    "concept_id": "hgnc:12012",
                },
                "elementGenomicStart": {
                    "id": "ga4gh:SL.Q8vkGp7_xR9vI0PQ7g1IvUUeQ4JlJG8l",
                    "digest": "Q8vkGp7_xR9vI0PQ7g1IvUUeQ4JlJG8l",
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
                    "end": 154192135,
                },
                "elementGenomicEnd": {
                    "id": "ga4gh:SL.Lnne0bSsgjzmNkKsNnXg98FeJSrDJuLb",
                    "digest": "Lnne0bSsgjzmNkKsNnXg98FeJSrDJuLb",
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
                },
            },
            {
                "type": "GeneElement",
                "gene": {
                    "type": "Gene",
                    "symbol": "ALK",
                    "concept_id": "hgnc:427",
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
                "strand": 1,
            },
            {"type": "MultiplePossibleGenesElement"},
        ],
        "regulatoryElement": {
            "type": "RegulatoryElement",
            "regulatoryClass": "promoter",
            "associatedGene": {
                "type": "Gene",
                "symbol": "BRAF",
                "concept_id": "hgnc:1097",
            },
        },
    }
