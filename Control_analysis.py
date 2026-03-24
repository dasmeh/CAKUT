import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# =========================
# USER INPUTS
# Edit these paths/parameters for your own dataset
# =========================
H5AD_FILE = "/path/to/input_data.h5ad"
CONTROL_FILE = "/path/to/control_table.csv"        # required columns: GENES, EXP, LENGTH
TARGET_GENE_FILE = "/path/to/target_genes.txt"     # one gene per line

CELLTYPE_COLUMN = "cell_type"                      # column in adata.obs
ASSAY_LAYER = None                                 # None = use adata.X
N_CONTROL_SETS = 1000
N_BINS_EXPR = 20
N_BINS_LEN = 20
MIN_GENES_PER_BIN = 5
RANDOM_SEED = 42
OUTPUT_PREFIX = "preferential_expression"

np.random.seed(RANDOM_SEED)


# =========================
# Helper functions
# =========================
def read_gene_list(path):
    genes = pd.read_csv(path, header=None)[0].astype(str).str.strip().tolist()
    genes = [g for g in genes if g]
    return list(dict.fromkeys(genes))


def get_expression_matrix(adata, layer=None):
    X = adata.layers[layer] if layer is not None else adata.X
    if sparse.issparse(X):
        return X
    return np.asarray(X)


def row_means_for_groups(X, groups):
    unique_groups = pd.Index(sorted(groups.unique()))
    result = []

    for g in unique_groups:
        idx = np.where(groups.values == g)[0]
        X_sub = X[idx, :]
        if sparse.issparse(X_sub):
            means = np.asarray(X_sub.mean(axis=0)).ravel()
        else:
            means = X_sub.mean(axis=0)
        result.append(means)

    mat = np.column_stack(result)
    return pd.DataFrame(mat, columns=unique_groups)


def overall_mean(X):
    if sparse.issparse(X):
        return np.asarray(X.mean(axis=0)).ravel()
    return X.mean(axis=0)


def make_bins(values, n_bins):
    return pd.qcut(values, q=n_bins, duplicates="drop")


def build_matched_pool(control_df, available_genes,
                       gene_col="GENES", expr_col="EXP", len_col="LENGTH",
                       n_bins_expr=20, n_bins_len=20):
    df = control_df.copy()
    df[gene_col] = df[gene_col].astype(str)

    df = df[df[gene_col].isin(available_genes)].dropna(subset=[expr_col, len_col]).copy()
    df = df.drop_duplicates(subset=gene_col).copy()

    df["expr_bin"] = make_bins(df[expr_col], n_bins_expr)
    df["len_bin"] = make_bins(df[len_col], n_bins_len)

    return df


def make_exact_matched_control_sets(target_genes, pool_df,
                                    gene_col="GENES",
                                    n_sets=1000,
                                    min_genes_per_bin=5):
    gene_to_bin = pool_df.set_index(gene_col)[["expr_bin", "len_bin"]].to_dict("index")
    valid_target_genes = [g for g in target_genes if g in gene_to_bin]

    if len(valid_target_genes) == 0:
        raise ValueError("None of the target genes were found in the matched pool.")

    grouped = {}
    for (eb, lb), subdf in pool_df.groupby(["expr_bin", "len_bin"], observed=True):
        grouped[(eb, lb)] = subdf[gene_col].tolist()

    control_sets = []

    for _ in range(n_sets):
        sampled = []
        for g in valid_target_genes:
            eb = gene_to_bin[g]["expr_bin"]
            lb = gene_to_bin[g]["len_bin"]
            candidates = grouped.get((eb, lb), [])

            candidates_excl = [x for x in candidates if x != g]

            if len(candidates_excl) >= min_genes_per_bin:
                chosen = np.random.choice(candidates_excl)
            elif len(candidates) > 0:
                chosen = np.random.choice(candidates)
            else:
                continue

            sampled.append(chosen)

        control_sets.append(sampled)

    return valid_target_genes, control_sets


def compute_preferential_score(mean_by_celltype_df, overall_means, gene_names, gene_set):
    gene_set = set(gene_set)
    idx = [i for i, g in enumerate(gene_names) if g in gene_set]

    if len(idx) == 0:
        raise ValueError("No genes from gene set found in adata.var_names")

    observed_by_ct = mean_by_celltype_df.iloc[idx, :].mean(axis=0)
    overall = overall_means[idx].mean()
    score = observed_by_ct - overall
    return score


def empirical_pvals(observed_scores, null_score_df):
    pvals = {}
    for ct in observed_scores.index:
        obs = observed_scores.loc[ct]
        null = null_score_df.loc[:, ct].values
        p = (1 + np.sum(null >= obs)) / (1 + len(null))
        pvals[ct] = p
    return pd.Series(pvals)


# =========================
# Load data
# =========================
print("Loading h5ad...")
adata = sc.read_h5ad(H5AD_FILE)

if CELLTYPE_COLUMN not in adata.obs.columns:
    raise ValueError(f"{CELLTYPE_COLUMN} not found in adata.obs")

print("Loading target genes...")
target_genes = read_gene_list(TARGET_GENE_FILE)

print("Loading control table...")
control_df = pd.read_csv(CONTROL_FILE)

required_cols = {"GENES", "EXP", "LENGTH"}
missing = required_cols - set(control_df.columns)
if missing:
    raise ValueError(f"Control file is missing columns: {missing}")

# Harmonize gene IDs
adata.var_names = adata.var_names.astype(str)
gene_names = adata.var_names.tolist()

# Optional: strip version numbers from Ensembl IDs like ENSG000001234.5
gene_names_clean = [g.split(".")[0] for g in gene_names]
adata.var["original_var_names"] = adata.var_names
adata.var_names = gene_names_clean
gene_names = adata.var_names.tolist()

# Intersect target genes with h5ad genes
target_genes_in_data = [g for g in target_genes if g in gene_names]

print(f"Target genes provided: {len(target_genes)}")
print(f"Target genes found in h5ad: {len(target_genes_in_data)}")

# =========================
# Expression summaries
# =========================
print("Computing group means...")
X = get_expression_matrix(adata, ASSAY_LAYER)
celltypes = adata.obs[CELLTYPE_COLUMN].astype(str)

mean_by_ct = row_means_for_groups(X, celltypes)
mean_by_ct.index = gene_names

overall = overall_mean(X)

# =========================
# Build matched gene pool
# =========================
print("Building matched pool...")
pool_df = build_matched_pool(
    control_df=control_df,
    available_genes=set(gene_names),
    gene_col="GENES",
    expr_col="EXP",
    len_col="LENGTH",
    n_bins_expr=N_BINS_EXPR,
    n_bins_len=N_BINS_LEN
)

target_genes_matched, control_sets = make_exact_matched_control_sets(
    target_genes=target_genes_in_data,
    pool_df=pool_df,
    gene_col="GENES",
    n_sets=N_CONTROL_SETS,
    min_genes_per_bin=MIN_GENES_PER_BIN
)

print(f"Target genes retained for matching: {len(target_genes_matched)}")
print(f"Generated {len(control_sets)} matched control sets")

# =========================
# Observed score
# =========================
print("Computing observed scores...")
observed_scores = compute_preferential_score(
    mean_by_celltype_df=mean_by_ct,
    overall_means=overall,
    gene_names=gene_names,
    gene_set=target_genes_matched
)

# =========================
# Null distribution
# =========================
print("Computing matched null scores...")
null_scores = []

for gs in control_sets:
    score = compute_preferential_score(
        mean_by_celltype_df=mean_by_ct,
        overall_means=overall,
        gene_names=gene_names,
        gene_set=gs
    )
    null_scores.append(score)

null_score_df = pd.DataFrame(null_scores)

# =========================
# Statistics
# =========================
print("Computing empirical p-values and z-scores...")
pvals = empirical_pvals(observed_scores, null_score_df)

null_mean = null_score_df.mean(axis=0)
null_std = null_score_df.std(axis=0, ddof=1).replace(0, np.nan)
z_scores = (observed_scores - null_mean) / null_std

results = pd.DataFrame({
    "cell_type": observed_scores.index,
    "observed_score": observed_scores.values,
    "null_mean": null_mean.reindex(observed_scores.index).values,
    "null_sd": null_std.reindex(observed_scores.index).values,
    "z_score": z_scores.reindex(observed_scores.index).values,
    "empirical_p": pvals.reindex(observed_scores.index).values,
    "n_target_genes_used": len(target_genes_matched),
}).sort_values(["empirical_p", "z_score"], ascending=[True, False])

# Benjamini-Hochberg FDR
results = results.sort_values("empirical_p").reset_index(drop=True)
m = len(results)
results["fdr_bh"] = results["empirical_p"] * m / np.arange(1, m + 1)
results["fdr_bh"] = np.minimum.accumulate(results["fdr_bh"][::-1])[::-1]
results = results.sort_values("z_score", ascending=False)

# =========================
# Save outputs
# =========================
results.to_csv(f"{OUTPUT_PREFIX}.celltype_results.tsv", sep="\t", index=False)
null_score_df.to_csv(f"{OUTPUT_PREFIX}.null_scores.tsv", sep="\t", index=False)
pd.Series(target_genes_matched, name="target_genes_used").to_csv(
    f"{OUTPUT_PREFIX}.target_genes_used.tsv",
    sep="\t",
    index=False
)

print("\nDone.")
print(results.head(20))