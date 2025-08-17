.. _install:

Installation
============

Prerequisites
-------------

* Python 3.11+
* A UNIX-like environment (e.g. MacOS, WSL, Ubuntu)
* A recent version of PostgreSQL (ideally at least 11+)
* A modern Java runtime (if using DynamoDB for the Gene Normalizer database)

Library installation
--------------------

Install ``FUSOR`` from `PyPI <https://pypi.org/project/fusor/>`_:

.. code-block:: shell

    pip install fusor

Data setup
----------

Universal Transcript Archive (UTA)
++++++++++++++++++++++++++++++++++

The `UTA <https://github.com/biocommons/uta>`_ is a dataset of genome-transcript aligned data supplied as a PostgreSQL database. Access in FUSOR is supplied by way of ``Cool-Seq-Tool``; see the `Cool-Seq-Tool UTA docs <https://coolseqtool.readthedocs.io/stable/install.html#set-up-uta>`_ for some opinionated setup instructions.

At runtime, UTA connection information can be relayed to FUSOR (by way of Cool-Seq-Tool) either as an initialization argument or via the environment variable ``UTA_DB_URL``. By default, it is set to ``postgresql://uta_admin:uta@localhost:5432/uta/uta_20210129b``. See the `Cool-Seq-Tool configuration docs <https://coolseqtool.readthedocs.io/stable/usage.html#environment-configuration>`_ for more info.

SeqRepo
+++++++

`SeqRepo <https://github.com/biocommons/biocommons.seqrepo>`_ is a controlled dataset of biological sequences. As with UTA, access in FUSOR is given via `Cool-Seq-Tool`, which provides `documentation <https://coolseqtool.readthedocs.io/stable/install.html#set-up-seqrepo>`_ on getting it set up.

At runtime, the file location of the SeqRepo instance directory can be defined (by way of Cool-Seq-Tool) either as an initialization argument or via the environment variable ``SEQREPO_ROOT_DIR``. By default, it's expected to be ``/usr/local/share/seqrepo/latest``. See the `Cool-Seq-Tool configuration docs <https://coolseqtool.readthedocs.io/stable/usage.html#environment-configuration>`_ for more info.

Gene Normalizer
+++++++++++++++

Finally, ``FUSOR`` uses the `Gene Normalizer <https://github.com/cancervariants/gene-normalization>`_ to ground gene terms. See the `Gene Normalizer documentation <https://gene-normalizer.readthedocs.io/stable/install.html>`_ for setup instructions.

Connection information for the normalizer database can be set using the environment variable ``GENE_NORM_DB_URL``. See the `Gene Normalizer docs <https://gene-normalizer.readthedocs.io/stable/reference/api/database/gene.database.database.html#gene.database.database.create_db>`_ for more information on connection configuration.
As a default, this connects to port 8000: ``http://localhost:8000``.

Check data availability
+++++++++++++++++++++++

Use the :py:meth:`fusor.tools.check_data_resources` method to verify that all data dependencies are available:

.. code-block:: pycon

   >>> from fusor.tools import check_data_resources
   >>> status = await check_data_resources()
   >>> assert all(status)  # passes if all resources can be acquired successfully
