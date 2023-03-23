import sys
import scanpy as sc
# import os.path

args = sys.argv
ad_out = args[1]
batch_column = str(args[2])
ad_path_list = args[3:]

# Load list/series/dict of anndata objects
ad_list = [sc.read_h5ad(h5ad_path) for h5ad_path in ad_path_list]
print(f"Loaded {len(ad_list)} anndata objects...")

# Get sample IDs for merging
# batch_categories = [os.path.splitext(os.path.basename(file_path))[0] for file_path in ad_path_list] # only works if file is /path/to/{sampleID}.h5ad
print(ad_list[0].obs_keys)
batch_categories = [ad.obs[batch_column][0] for ad in ad_list]
print(batch_categories)

# Merge
merged_ad = ad_list[0].concatenate(
    ad_list[1:],
    uns_merge="same",
    index_unique="-",
    batch_categories=batch_categories
)
merged_ad.obs_names_make_unique()
print(f"Merged object has {merged_ad.shape[0]} cells and {merged_ad.shape[1]} genes")

# Write merged anndata object
merged_ad.write(
    ad_out
)