.. data_model:

Data Model
==========

FUSOR's model classes implement the objects described in the `Minimum Information Model <https://fusions.cancervariants.org/en/latest/information_model.html>`_ of the `VICC Gene Fusion Specification <https://fusions.cancervariants.org/en/latest/index.html>`_. while we aim to promptly adopt new features and modifications made to the Gene Fusion Specification, as of June 2025, it does not yet have discrete releases. Therefore, FUSOR utilizes the latest version of the Specification as described in its documentation. Once the Specification is stable, we will align specific FUSOR releases with discrete VICC Fusion Specification versions.

FUSOR is an *implementation* of the Minimum Information Model, meaning it introduces some additional opinionated structure on top of what the specification requires. In particular, it employs a number of genomic description and location classes as defined in the `GA4GH Variation Representation Specification (VRS) <https://vrs.ga4gh.org/en/stable/>`_ and implemented in `VRS-Python <https://github.com/ga4gh/vrs-python>`_ version |vrs_python_version|.
