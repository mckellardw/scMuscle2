import sys
from scanpy import read_10x_mtx

args = sys.argv
# print(args)
mat_in = args[1]
ad_out = args[2]

adata = read_10x_mtx(
    path=mat_in,
    var_names="gene_symbols",
    make_unique=True,
    cache=False
)

adata.write(ad_out)
