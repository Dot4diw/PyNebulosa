"""
Example usage of PyNebulosa for gene expression density plotting
"""

import numpy as np
import matplotlib.pyplot as plt
import pynebulosa as pyn
import scanpy as sc

adata = sc.read("./datasets/pancreas.h5ad")

pyn.plot_density(adata, 
                 basis="umap",
                 features=["Spp1","Rps19"],
                 joint=True, 
                 raster= False, 
                 legend_shrink=1.0,
                 pal="plasma",
                 figsize = (4.5,3.5))
