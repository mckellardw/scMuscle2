import sys
import scanpy as sc

#TODO- change this from positional arguments...
args = sys.argv
ad_in = args[1]
ad_out = args[2]
min_genes = args[3]
min_counts = args[4]
max_mito = args[5]

adata = sc.read_h5ad(ad_in)

sc.pp.filter_cells(
    adata,
    min_genes=min_genes
)
sc.pp.filter_cells(
    adata,
    min_counts=min_counts
)

# annotate the group of mitochondrial genes as 'mito'
adata.var['mito'] = adata.var_names.str.startswith('MT-')

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=['mito'],
    percent_top=None,
    log1p=False,
    inplace=True
)

adata = adata[adata.obs.pct_counts_mito < max_mito, :]

# Doublet prediction/filtering
sc.external.pp.scrublet(
    adata
)
# sc.external.pl.scrublet_score_distribution(
#     adata,
#     figsize =[6,2.25]
# )
cutoff_threshold = 42 #TODO- algorithmically compute this...
adata = adata[adata.obs["doublet_score"] < cutoff_threshold,]


adata.write(ad_out)
