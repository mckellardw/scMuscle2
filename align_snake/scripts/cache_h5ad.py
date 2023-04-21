import sys
from scanpy import read_10x_mtx

args = sys.argv
mat_in = args[1]
ad_out = args[2]
var_names = args[3]

# Read in counts
adata = read_10x_mtx(
    path=mat_in,
    var_names=var_names,
    make_unique=True,
    cache=False
)

# Save raw count matrix
adata.raw = adata

adata.write(
    ad_out,
    compression='gzip'
)
