"""Provide miscellaneous tools for fusion modeling."""

import logging
from collections import namedtuple
from typing import Annotated

from biocommons.seqrepo.seqrepo import SeqRepo
from cool_seq_tool.app import CoolSeqTool
from cool_seq_tool.resources.status import check_status as check_cst_status
from gene.database import AbstractDatabase as GeneDatabase
from gene.database import create_db
from gene.schemas import CURIE_REGEX
from pydantic import StringConstraints, ValidationError

from fusor.exceptions import IDTranslationException

_logger = logging.getLogger(__name__)


def translate_identifier(
    seqrepo: SeqRepo, ac: str, target_namespace: str = "ga4gh"
) -> Annotated[str, StringConstraints(pattern=CURIE_REGEX)]:
    """Return ``target_namespace`` identifier for accession provided.

    :param ac: Identifier accession
    :param target_namespace: The namespace of identifiers to return.
        Default is ``ga4gh``
    :return: Identifier for ``target_namespace``
    :raise: IDTranslationException if unable to perform desired translation
    """
    try:
        target_ids = seqrepo.translate_identifier(
            ac, target_namespaces=target_namespace
        )
    except KeyError as e:
        _logger.warning("Unable to get translated identifier: %s", e)
        raise IDTranslationException from e

    if not target_ids:
        raise IDTranslationException
    return target_ids[0]


FusorDataResourceStatus = namedtuple(
    "FusorDataResourceStatus", ("cool_seq_tool", "gene_normalizer")
)


async def check_data_resources(
    gene_database: GeneDatabase | None = None,
    cool_seq_tool: CoolSeqTool | None = None,
) -> FusorDataResourceStatus:
    """Perform basic status checks on known data requirements.

    Mirroring the input structure of the :py:class:`fusor.fusor.FUSOR` class, existing
    instances of the Gene Normalizer database and Cool-Seq-Tool can be passed as
    arguments. Otherwise, resource construction is attempted in the same manner as it
    would be with the FUSOR class, relying on environment variables and defaults.

    >>> from fusor.tools import check_data_resources
    >>> status = await check_data_resources()
    >>> assert all(status)  # passes if all resources can be acquired successfully

    The return object is a broad description of resource availability, grouped by
    library. For a more granular description to support debugging, all failures are
    logged as ``logging.ERROR`` by respective upstream libraries.

    :param gene_database: gene normalizer DB instance
    :param cool_seq_tool: Cool-Seq-Tool instance
    :return: namedtuple describing whether Cool-Seq-Tool and Gene Normalizer resources
        are all available
    """
    if cool_seq_tool is None:
        cool_seq_tool = CoolSeqTool()
    cst_status = await check_cst_status()

    if gene_database is None:
        gene_database = create_db()

    gene_status = False
    try:
        if not gene_database.check_schema_initialized():
            _logger.error("Health check failed: gene DB schema uninitialized")
        else:
            if not gene_database.check_tables_populated():
                _logger.error("Health check failed: gene DB is incompletely populated")
            else:
                gene_status = True
    except Exception as e:
        _logger.error("Encountered error while creating gene DB: %s", e)
    return FusorDataResourceStatus(
        cool_seq_tool=all(cst_status), gene_normalizer=gene_status
    )


def get_error_message(e: ValidationError) -> str:
    """Get all error messages from a pydantic ValidationError

    :param e: the ValidationError to get the messages from
    :return: string containing all of the extracted error messages, separated by newlines or the string
    representation of the exception if 'msg' field is not present
    """
    if e.errors():
        return "\n".join(str(error["msg"]) for error in e.errors() if "msg" in error)
    return str(e)
