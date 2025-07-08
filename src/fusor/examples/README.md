# Examples

* `bcr_abl1.json`: Example BCR-ABL1 categorical fusion drawn from [COSF1780](https://cancer.sanger.ac.uk/cosmic/fusion/summary?id=1780). Demonstrates structure of junction components, a linker sequence segment, critical functional domains, and reading frame preservation. Represented in nomenclature as `NM_004327.3(BCR):e.2+182::ACTAAAGCG::NM_005157.5(ABL1):e.2-173`.
* `bcr_abl1_expanded.json`: Equivalent fusion to the above, but with expanded descriptions of genes, locations, and sequences provided by SeqRepo and the VICC Gene Normalizer.
* `alk.json`: Example of an ALK fusion, demonstrating use of a categorical "multiple possible gene" component, retrieved from a human-curated database like [CIViC](https://civicdb.org/variants/499/summary). Represented in nomenclature as `v::ALK(hgnc:427)`.
* `ewsr1.json`: An EWSR1 assayed fusion, demonstrating an assay description object and use of the "unknown gene" partner. Represented in nomenclature as `EWSR1(hgnc:3508)::?`.
* `tpm3_itd.json`: Example of an Internal Tandem Duplication involving exons 1-8 of TPM3. There is no nomenclature string as ITDs are not supported in the VICC Gene Fusion Specification.
* `tpm3_ntrk1.json`: Example TPM3-NTRK1 assayed fusion drawn from previous VICC Fusion Curation draft material. Represented in nomenclature as `NM_152263.3(TPM3):e.1_8::NM_002529.3(NTRK1):e.10_22`.
* `tpm3_pdgfrb.json`: Example TPM3-PDGFRB assayed fusion identified via RT-PCR. Represented in nomenclature as `NM_152263.3(TPM3):e.8::NM_002609.3(PDGFRB):e.11_22`.
* `igh_myc.json`: Example of an enhancer-driven IGH-MYC categorical fusion. Represented in nomenclature as `reg_e_EH38E3121735@IGH(hgnc:5477)::MYC(hgnc:7553)`.
