
def knee_plot(
    ADATA,
    x_lim=[0, 20000],
    line_width=2,
    line_color="b"
    verbose=False
):
    import matplotlib.pyplot as plt
    expected_num_cells = 10000

    for i in range(0,meta.shape[0]):
        knee = np.sort((np.array(ADATA.X.sum(axis=1))).flatten())[::-1]

        fig, ax = plt.subplots(figsize=(10, 7))

        ax.loglog(
            knee,
            range(len(knee)),
            linewidth=line_width,
            color=line_color
        )
    #     ax.axvline(x=knee[expected_num_cells], linewidth=3, color="k")
    #     ax.axhline(y=expected_num_cells, linewidth=3, color="k")

        ax.set_xlabel("UMI Counts")
        ax.set_ylabel("Set of Barcodes")
        ax.set_title(ADATA.obs["sample"][0])

        plt.xlim(x_lim)
        plt.grid(True, which="both")
        plt.show()
