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
FUSOR relies on `Seqrepo <https://github.com/biocommons/biocommons.seqrepo>`_, which you must download yourself.

FUSOR uses Seqrepo to retrieve sequences at given positions on a transcript.

From the *root* directory:

.. code-block:: shell

    pip install seqrepo
    sudo mkdir /usr/local/share/seqrepo
    sudo chown $USER /usr/local/share/seqrepo
    seqrepo pull -i 2024-12-20/  # Replace with latest version using `seqrepo list-remote-instances` if outdated

If you get an error similar to the one below:

.. code-block:: shell

    PermissionError: [Error 13] Permission denied: '/usr/local/share/seqrepo/2024-12-20/._fkuefgd' -> '/usr/local/share/seqrepo/2024-12-20/'

You will want to do the following:
*(_Might not be .\_fkuefgd, so replace with your error message path_)*

.. code-block:: shell

    sudo mv /usr/local/share/seqrepo/2024-12-20._fkuefgd /usr/local/share/seqrepo/2024-12-20
    exit

Use the ``SEQREPO_ROOT_DIR`` environment variable to set the path of an already existing SeqRepo directory. The default is ``/usr/local/share/seqrepo/latest``.

Gene Normalizer
+++++++++++++++

Finally, ``FUSOR`` uses the `Gene Normalizer <https://github.com/cancervariants/gene-normalization>`_ to ground gene terms. See the `Gene Normalizer documentation <https://gene-normalizer.readthedocs.io/stable/install.html>`_ for setup instructions.

Connection information for the normalizer database can be set using the environment variable ``GENE_NORM_DB_URL``. See the `Gene Normalizer docs <https://gene-normalizer.readthedocs.io/stable/reference/api/database/gene.database.database.html#gene.database.database.create_db>`_ for more information on connection configuration.
As a default, this connects to port 8000: ``http://localhost:8000``.

Docker
++++++

FUSOR's dependencies can be installed using a Docker container.

.. important::

   This section assumes you have a local
   `SeqRepo <https://github.com/biocommons/biocommons.seqrepo>`_
   installed at ``/usr/local/share/seqrepo/2024-12-20``.

   If you're using Docker Desktop, you must go to
   **Settings → Resources → File sharing** and add
   ``/usr/local/share/seqrepo`` under the *Virtual file shares*
   section. Otherwise, you will get the following error::

      OSError: Unable to open SeqRepo directory /usr/local/share/seqrepo/2024-12-20

To build, (re)create, and start containers:

.. code-block:: shell

   docker volume create uta_vol
   docker compose up

.. tip::

   If you want a clean slate, run ``docker compose down -v`` to remove
   containers and volumes, then run
   ``docker compose up --build`` to rebuild and start fresh containers.

In Docker Desktop, you should see the following for a successful setup:

.. figure:: ../../docker-desktop-container.png
   :alt: Docker Desktop Container
   :align: center

.. note::

   `python-dotenv <https://pypi.org/project/python-dotenv/>`_ can be used
   to load environment variables needed for analysis notebooks in the
   ``notebooks`` directory. Environment variables can be found in
   ``.env.shared``.

Check data availability
+++++++++++++++++++++++

Use the :py:meth:`fusor.tools.check_data_resources` method to verify that all data dependencies are available:

.. code-block:: pycon

   >>> from fusor.tools import check_data_resources
   >>> status = await check_data_resources()
   >>> assert all(status)  # passes if all resources can be acquired successfully
