# Utility functions for use with scanpy

# Calculate the number of PCs that contain some proportion (default is 95%) of the variance
def npcs(
  ADATA,
  var_perc=0.95,
  reduction="pca"
):
    from numpy import sum, var
    get_var = lambda i: var(ADATA.obsm[reduction][:,i])

    if ADATA.obsm[reduction] is None:
        print(f"Reduction '{reduction}', not found!")
        return None
    else:
        var_tmp = [get_var(i) for i in list(range(0,ADATA.obsm[reduction].shape[1]))]
        var_cut = var_perc * sum(var_tmp)
        n_pcs = 0
        var_sum = 0
        while var_sum<var_cut and n_pcs<ADATA.obsm[reduction].shape[1]-1:
            var_sum = var_sum + var_tmp[n_pcs]
            n_pcs = n_pcs + 1

        return(n_pcs)


# Re-order a dimensions of a reduction by decreasing % variance
def reorder_reduction(
    ADATA,
    reduction="pca",
    verbose=False
):
    from numpy import var, argsort

    if reduction in ADATA.obsm:
        get_var = lambda i: var(ADATA.obsm[reduction][:,i])
        var_tmp = [get_var(i) for i in list(range(0,ADATA.obsm[reduction].shape[1]))]
        if verbose:
            print("Reduction variance by dimension:")
            print(var_tmp)

        pc_order = argsort(var_tmp)[::-1]
        ADATA.obsm[reduction] = ADATA.obsm[reduction][:,pc_order]
    else:
        print(f"The reduction '{reduction}' was not found...")
