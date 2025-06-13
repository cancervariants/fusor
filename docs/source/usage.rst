.. _usage:

Usage
=====

Constructing fusions
--------------------

:py:class:`fusor.fusor.FUSOR` is the core class of the FUSOR package, enabling easy construction of new fusions given relevant parameters:

.. code-block:: pycon

   >>> from fusor import FUSOR
   >>> f = FUSOR()
   >>> fusion = f.fusion(
   ...     structure=[
   ...         {
   ...             "type": "GeneElement",
   ...             "gene": {
                        "primaryCoding": {
                              "id": "hgnc:3508",
                              "code": "HGNC:3508",
                              "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/"
                        },
                        "conceptType": "Gene",
                        "name": "EWSR1"
                  }
   ...         },
   ...         {
   ...           "type": "UnknownGeneElement"
   ...         }
   ...     ],
   ...     assay={
   ...         "type": "Assay",
   ...         "methodUri": "pmid:33576979",
   ...         "assayId": "obi:OBI_0003094",
   ...         "assayName": "fluorescence in-situ hybridization assay",
   ...         "fusionDetection": "inferred",
   ...     },
   ... )
   >>> fusion.type
   <FUSORTypes.ASSAYED_FUSION: 'AssayedFusion'>

As seen in the example above, :py:meth:`fusor.fusor.FUSOR.fusion()` can infer fusion context (i.e. whether it is `categorical or assayed <https://fusions.cancervariants.org/en/latest/terminology.html#gene-fusion-contexts>`_), assuming the provided arguments are sufficient for Pydantic to discern which kind of fusion it would have to be. In cases where the type is inherently ambiguous, an explicit type parameter can be passed, or the fusion's ``type`` property can be included:

.. code-block:: pycon

   >>> from fusor import FUSOR
   >>> f = FUSOR()
   >>> structure = [
           {"gene": {"primaryCoding": {"id": "hgnc:3508", "code": "HGNC:3508", "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/"}, "conceptType": "Gene", "name": "EWSR1"},  "type": "GeneElement"},
           {"gene": {"primaryCoding": {"id": "hgnc:3446", "code": "HGNC:3446", "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/"}, "conceptType": "Gene", "name": "ERG"},  "type": "GeneElement"}
   ... ]
   >>> fusion = f.fusion(**{"structure": structure, "type": "CategoricalFusion"})
   >>> fusion.type
   <FUSORTypes.Categorical_FUSION: 'CategoricalFusion'>
   >>> from fusor.models import FusionType
   >>> fusion = f.fusion(fusion_type=FusionType.Categorical_FUSION, **{"structure": structure})
   >>> fusion.type
   <FUSORTypes.CATEGORICAL_FUSION: 'CategoricalFusion'>

Validating fusions
------------------

``fusor.models`` defines classes for fusions, as well as the associated subcomponents described in the `fusions minimal information model <https://fusions.cancervariants.org/en/latest/information_model.html>`_. These are implemented as `Pydantic <https://docs.pydantic.dev/latest/>`_ classes, enabling runtime type validation.

.. code-block:: pycon

   >>> from fusor.models import GeneElement, UnknownGeneElement, AssayedFusion
   >>> gene = GeneElement(gene={"primaryCoding": {"id": "hgnc:1097", "code": "HGNC:1097", "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/"}, "conceptType": "Gene", "name": "BRAF"})
   >>> fusion = AssayedFusion(structure=[UnknownGeneElement(), gene])
   >>> fusion
   AssayedFusion(type=<FUSORTypes.ASSAYED_FUSION: 'AssayedFusion'>, regulatoryElement=None, structure=[UnknownGeneElement(type=<FUSORTypes.UNKNOWN_GENE_ELEMENT: 'UnknownGeneElement'>), GeneElement(type=<FUSORTypes.GENE_ELEMENT: 'GeneElement'>, gene=Gene(primaryCoding=(id="hgnc:1097", "code": "HGNC:1097", "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/"}, "conceptType": "Gene", "name": "BRAF")))], readingFramePreserved=None, causativeEvent=None, assay=None)

In this example, a fusion is constructed with only one structural element, even though fusions are defined as `"the joining of two or more genes" <https://fusions.cancervariants.org/en/latest/terminology.html#gene-fusions>`_:

.. code-block:: pycon

   >>> from fusor.models import AssayedFusion
   >>> AssayedFusion(**{"structure": [{"gene": {"primaryCoding": {"id": "hgnc:3508", "code": "HGNC:3508", "system": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/"}, "conceptType": "Gene", "name": "EWSR1"},  "type": "GeneElement"}]})
   Traceback (most recent call last):
     File "<stdin>", line 1, in <module>
     File "/Users/jss009/code/fusor/.venv/lib/python3.12/site-packages/pydantic/main.py", line 159, in __init__
       __pydantic_self__.__pydantic_validator__.validate_python(data, self_instance=__pydantic_self__)
   pydantic_core._pydantic_core.ValidationError: 1 validation error for AssayedFusion
     Value error, Fusions must contain >= 2 structural elements, or >=1 structural element and a regulatory element [type=value_error, input_value={'structure': [{'type': '...', 'name': 'EWSR1'}}]}, input_type=dict]
       For further information visit https://errors.pydantic.dev/2.1/v/value_error

Example fusions
---------------

``fusor.examples`` contains pre-defined fusion objects intended to illustrate various aspects of the information model and nomenclature. ``fusor.examples.alk`` represents the category of fusions between the ALK gene and any other partner, where the `protein kinase, ATP binding site domain <https://www.ebi.ac.uk/interpro/entry/InterPro/IPR017441/>`_ is preserved:

.. code-block:: pycon

   >>> from fusor import examples
   >>> examples.alk.type
   <FUSORTypes.CATEGORICAL_FUSION: 'CategoricalFusion'>
   >>> examples.alk.structure[0]
   MultiplePossibleGenesElement(type=<FUSORTypes.MULTIPLE_POSSIBLE_GENES_ELEMENT: 'MultiplePossibleGenesElement'>)
   >>> examples.alk.structure[1]
   GeneElement(type=<FUSORTypes.GENE_ELEMENT: 'GeneElement'>, gene=Gene(primaryCoding=(id="hgnc:427", code="HGNC:427", system="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/"), conceptType="Gene", name="ALK"))
   >>> examples.alk.criticalFunctionalDomains[0].status
   <DomainStatus.PRESERVED: 'preserved'>
   >>> examples.alk.criticalFunctionalDomains[0].id
   'interpro:IPR017441'

Provided examples include:

* ``examples.bcr_abl1``: Example BCR-ABL1 categorical fusion drawn from `COSF1780 <https://cancer.sanger.ac.uk/cosmic/fusion/summary?id=1780>`_. Demonstrates structure of junction components, a linker sequence segment, critical functional domains, and reading frame preservation. Represented in nomenclature as ``NM_004327.3(BCR):e.2+182::ACTAAAGCG::NM_005157.5(ABL1):e.2-173``.
* ``examples.bcr_abl1_expanded``: Equivalent fusion to the above, but with expanded descriptions of genes, locations, and sequences provided by SeqRepo and the VICC Gene Normalizer.
* ``examples.alk``: Example of an ALK fusion, demonstrating use of a categorical "multiple possible gene" component, retrieved from a human-curated database like `CIViC <https://civicdb.org/variants/499/summary>`_. Represented in nomenclature as ``v::ALK(hgnc:427)``.
* ``examples.ewsr1``: An EWSR1 assayed fusion, demonstrating an assay description object and use of the "unknown gene" partner. Represented in nomenclature as ``EWSR1(hgnc:3508)::?``.
* ``examples.tpm3_ntrk1``: Example TPM3-NTRK1 assayed fusion drawn from previous VICC Fusion Curation draft material. Represented in nomenclature as ``NM_152263.3(TPM3):e.1_8::NM_002529.3(NTRK1):e.10_22``.
* ``examples.tpm3_pdgfrb``: Example TPM3-PDGFRB assayed fusion identified via RT-PCR. Represented in nomenclature as ``NM_152263.3(TPM3):e.8::NM_002609.3(PDGFRB):e.11_22``.
* ``examples.igh_myc``: Example of an enhancer-driven IGH-MYC categorical fusion. Represented in nomenclature as ``reg_e_EH38E3121735@IGH(hgnc:5477)::MYC(hgnc:7553)``.


Generating nomenclature
-----------------------

The core :py:class:`fusor.fusor.FUSOR` class can generate nomenclature for a fusion instance in line with the `VICC fusion nomenclature <https://fusions.cancervariants.org/en/latest/nomenclature.html>`_:

.. code-block:: pycon

   >>> from fusor import FUSOR, examples
   >>> f = FUSOR()
   >>> f.generate_nomenclature(examples.alk)
   'v::ALK(hgnc:427)'
