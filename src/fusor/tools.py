"""Provide miscellaneous tools for fusion modeling."""

import logging
from collections import namedtuple

from cool_seq_tool.resources.status import check_status as check_cst_status
from gene.database import AbstractDatabase as GeneDatabase
from gene.database import create_db
from pydantic import ValidationError

_logger = logging.getLogger(__name__)


FusorDataResourceStatus = namedtuple(
    "FusorDataResourceStatus", ("cool_seq_tool", "gene_normalizer")
)


async def check_data_resources(
    gene_database: GeneDatabase | None = None,
) -> FusorDataResourceStatus:
    """Perform basic status checks on known data requirements.

    Mirroring the input structure of the :py:class:`fusor.fusor.FUSOR` class, existing
    instances of the Gene Normalizer database can be passed as
    arguments. Otherwise, resource construction is attempted in the same manner as it
    would be with the FUSOR class, relying on environment variables and defaults.

    >>> from fusor.tools import check_data_resources
    >>> status = await check_data_resources()
    >>> assert all(status)  # passes if all resources can be acquired successfully

    The return object is a broad description of resource availability, grouped by
    library. For a more granular description to support debugging, all failures are
    logged as ``logging.ERROR`` by respective upstream libraries.

    :param gene_database: gene normalizer DB instance
    :return: namedtuple describing whether Cool-Seq-Tool and Gene Normalizer resources
        are all available
    """
    cst_status = await check_cst_status()

    gene_status = False
    try:
        if gene_database is None:
            gene_database = create_db()
        if not gene_database.check_schema_initialized():
            _logger.error("Health check failed: gene DB schema uninitialized")
        else:
            if not gene_database.check_tables_populated():
                _logger.error("Health check failed: gene DB is incompletely populated")
            else:
                gene_status = True
    except Exception:
        _logger.exception(
            "Encountered error while creating gene DB during resource check"
        )
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
