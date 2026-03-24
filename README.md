# Single-cell preferential expression of CAKUT-associated genes

**Congenital anomalies of the kidney and urinary tract (CAKUT)** comprise a heterogeneous spectrum of developmental disorders, yet it remains unclear whether their diverse genetic causes converge on shared cellular programs.

Here, we analyzed the expression of **91 CAKUT-associated genes** across single-cell RNA-seq datasets from **human fetal and adult kidney, ureter, and bladder**, complemented by **early embryonic gene expression data**.

During kidney development, canonical CAKUT genes such as **EYA1, SIX1, PAX2, and FOXC1** showed strong preferential expression in **mesangial** and **mesonephric nephron tubule epithelial cells**, highlighting early roles in **ureteric bud induction** and **branching morphogenesis**.

Temporal analyses revealed **two distinct expression trajectories**:

- **Early-peak genes**, including *EYA1, SIX1, SIX2, PAX2,* and *TGFA*, which peak during early nephrogenesis and decline thereafter  
- **Late-rise genes**, such as *MUC1*, which increase toward adult stages  

Together, these findings reveal a **conserved stromal–mesenchymal gene signature** underlying CAKUT pathogenesis and suggest disruption of **mesodermal lineage programs** that precede ureteric bud development and organ patterning.

---

## Methodology: scDRS-based cell-type prioritization

Cell-type–specific enrichment of CAKUT-associated genes was quantified using the **single-cell Disease Relevance Score (scDRS)** framework, which integrates disease gene sets with single-cell transcriptomic profiles while controlling for gene expression confounders.

**scDRS paper**  
https://www.nature.com/articles/s41588-022-01167-z

**scDRS GitHub repository**  
https://github.com/martinjzhang/scDRS

---

## Data availability

Processed single-cell datasets used in this study are provided as `.h5ad` files:

**Download processed h5ad files**  
https://www.dropbox.com/scl/fo/42v8bi6hk2qjozl5n5zr1/AMIiTv1r0hGRzkGVPVxdRok?rlkey=bnaeilegbmfzqq6ful5idmgey&dl=0

Datasets include:
- Human fetal kidney
- Adult kidney
- Ureter
- Bladder
- Early embryonic reference data

---

## Interactive data exploration (CELLxGENE)

The primary datasets used in this systematic analysis are publicly available and were obtained from CELLxGENE Explorer, where they can be interactively visualized and explored.

### Human fetal kidney
Developmental single-cell atlas capturing nephrogenesis, stromal and mesenchymal lineages.

https://cellxgene.cziscience.com/collections/bcb61471-2a44-4d00-a0af-ff085512674c

---

### Adult human kidney
Reference single-cell atlas of mature nephron and interstitial cell populations.

https://cellxgene.cziscience.com/collections/981b0c87-76d0-4edf-b201-58c52b5cd2f0

---

### Human ureter
Single-cell transcriptomic atlas of ureteral epithelial and stromal compartments.

https://cellxgene.cziscience.com/collections/0e54d4de-44f0-4d50-8649-b5c2bbe8f5d1

---

### Human bladder
Single-cell atlas of bladder urothelial, stromal, and immune populations.

https://cellxgene.cziscience.com/collections/af286675-7167-4816-a907-c72d41b73d37

---

## Reproducibility

- Analyses performed using Python and Scanpy-compatible workflows  
- CAKUT gene set curated from published literature  
- scDRS parameters and preprocessing steps documented in analysis scripts  

For questions or reuse, please open an issue.

---

## Citation

```bibtex
@misc{cakut_scRNAseq_2025,
  author = {Dasmeh, Pouria},
  title = {Single-cell preferential expression of CAKUT-associated genes},
  year = {2025},
  url = {https://github.com/<your-username>/<repo-name>}
}
