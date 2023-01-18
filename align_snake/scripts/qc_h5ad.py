import sys
import scanpy as sc

args = sys.argv
ad_in = args[1]
ad_out = args[2]
min_genes = int(args[3])
min_counts = int(args[4])
max_mito = int(args[5])

adata = sc.read_h5ad(
    ad_in
)

sc.pp.filter_cells(
    adata,
    min_genes=min_genes
)
sc.pp.filter_cells(
    adata,
    min_counts=min_counts
)

#preprocess according to Seurat recipe
# sc.pp.recipe_seurat(
#     adata
# )

# annotate the group of mitochondrial genes as 'mito'
#TODO- account for different species gene names...
adata.var['mito'] = adata.var_names.str.startswith('mt-')

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=['mito'],
    percent_top=None,
    log1p=False,
    inplace=True
)

# adata = adata[adata.obs.pct_counts_mito < max_mito, :]

# Doublet prediction/filtering
sc.external.pp.scrublet(
    adata
)
# sc.external.pl.scrublet_score_distribution(
#     adata,
#     figsize =[6,2.25]
# )

adata.write(ad_out)
