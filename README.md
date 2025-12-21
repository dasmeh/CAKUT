# Single-cell preferential expression of CAKUT-associated genes

**Congenital anomalies of the kidney and urinary tract (CAKUT)** comprise a heterogeneous spectrum of developmental disorders, yet it remains unclear whether their diverse genetic causes converge on shared cellular programs.

Here, we analyzed the expression of **91 CAKUT-associated genes** across single-cell RNA-seq datasets from **human fetal and adult kidney, ureter, and bladder**, complemented by **early embryonic gene expression data**.

During kidney development, canonical CAKUT genes such as **EYA1, SIX1, PAX2, and FOXC1** showed strong preferential expression in **mesangial** and **mesonephric nephron tubule epithelial cells**, highlighting early roles in **ureteric bud induction** and **branching morphogenesis**.

Temporal analyses further revealed **two distinct expression trajectories**:

- **Early-peak genes**, including *EYA1, SIX1, SIX2, PAX2,* and *TGFA*, which reached maximal expression during early nephrogenesis and declined thereafter  
- **Late-rise genes**, such as *MUC1*, which increased toward adult stages  

Together, these findings reveal a **conserved stromal–mesenchymal gene signature** underlying CAKUT pathogenesis and suggest that genetic perturbations may disrupt **mesodermal lineage programs** that prefigure ureteric bud development and organ patterning.

---

## Methodology: scDRS-based cell-type prioritization

Cell-type–specific enrichment of CAKUT-associated genes was quantified using the **single-cell Disease Relevance Score (scDRS)** framework, which integrates gene sets with single-cell transcriptomic profiles to identify disease-relevant cell populations while controlling for gene expression confounders.

📄 **scDRS paper**  
> Cao et al., *Cell*, 2023  
> https://www.cell.com/cell/fulltext/S0092-8674(22)01534-4  

📦 **scDRS GitHub repository**  
https://github.com/martinjzhang/scDRS

---

## Data availability

All processed single-cell datasets used in this analysis are provided as `.h5ad` files and are publicly available:

🔗 **Download single-cell h5ad files**  
https://www.dropbox.com/scl/fo/42v8bi6hk2qjozl5n5zr1/AMIiTv1r0hGRzkGVPVxdRok?rlkey=bnaeilegbmfzqq6ful5idmgey&dl=0

These datasets include:
- Human fetal kidney
- Adult kidney
- Ureter
- Bladder
- Early embryonic reference data

---

## Reproducibility

- Analysis performed using Python / Scanpy-compatible workflows  
- Gene prioritization based on curated CAKUT gene lists from the literature  
- scDRS parameters and preprocessing steps are documented in the analysis scripts  

For questions or reuse, please open an issue or contact the authors.

---

**Citation** (if using this resource):
```bibtex
@misc{cakut_scRNAseq_2025,
  author = {Dasmeh, Pouria},
  title = {Single-cell preferential expression of CAKUT-associated genes},
  year = {2025},
  url = {https://github.com/<your-username>/<repo-name>}
}
