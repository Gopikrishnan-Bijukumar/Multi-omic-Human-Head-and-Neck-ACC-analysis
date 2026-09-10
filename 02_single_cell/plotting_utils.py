import matplotlib.pyplot as plt
import scanpy as sc

def plot_umap(adata, color, s1, s2, size, out):
    fig, ax = plt.subplots(figsize = (int(s1), int(s2)))
    a = sc.pl.umap(
        adata,
        color = [str(color)],
        frameon = False,
        size = int(size),
        ax = ax,
        show = out
    )
    return a