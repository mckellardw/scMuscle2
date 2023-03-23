import sys
import scanpy as sc
import pandas as pd

args = sys.argv
ad_in = args[1]  # anndata file w/ out metadata
ad_out = args[2] # anndata file w/ metadata
meta = args[3]   # metadata (dataframe) to add

# Read in anndata
adata = sc.read_h5ad(ad_in)

# Open metadata .csv
meta_df = pd.read_csv(meta)

# Add metadata
for j in range(0,meta_df.shape[1]): 
    adata.obs[meta_df.columns[j]] = meta_df.iloc[0,j]

# Write new anndata
adata.write(ad_out)
