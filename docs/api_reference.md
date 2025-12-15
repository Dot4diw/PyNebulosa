# PyNebulosa API Reference

PyNebulosa is a Python implementation of the R Nebulosa package for visualizing gene expression data in single-cell RNA sequencing datasets. It provides kernel density estimation-based visualization methods that help recover signals from sparse data.

## Main Functions

### `plot_density`

```python
plot_density(adata, features, ax=None, basis=None, layer=None, dims=[1, 2], joint=False, only_joint_plot=False, method='wkde', adjust=1, shape='o', size=1, pal="viridis", raster=True, figsize=None, show_legend_title=True, legend_title="Density", legend_shrink=0.8, **kwargs)
```

Main function for plotting gene expression density maps.

#### Parameters

**Required Parameters:**
- `adata` (anndata.AnnData): AnnData object with cell embeddings and gene expression data
- `features` (str or list of str): Features (e.g. genes) to visualize

**Optional Parameters:**
- `ax` (matplotlib.axes.Axes, optional): Axes object to plot on
- `basis` (str, optional): Name of the embedding to visualize (e.g., 'umap', 'pca'). If not provided, uses the last embedding in .obsm
- `layer` (str, optional): Layer to use for expression data. If None, uses .X
- `dims` (list of int): Vector of length 2 specifying the dimensions to be plotted (default: [1, 2])
- `joint` (bool): Return joint density plot? By default False
- `only_joint_plot` (bool): When True and joint=True, only plot the joint density map (ignore individual gene plots). By default False
- `method` (str): Kernel density estimation method ('wkde'), (default: 'wkde')
- `adjust` (float): Numeric value to adjust to bandwidth (default: 1)
- `shape` (str or int): Shape of the points (default: 'o')
- `size` (float): Size of the points (default: 1)
- `pal` (str): Color palette name (default: "viridis")
- `raster` (bool): Rasterize plot (default: True)
- `figsize` (tuple of float, optional): Figure size as (width, height)
- `show_legend_title` (bool): Whether to show legend title. By default True
- `legend_title` (str): String used as legend title. Used when show_legend_title=True (default: "Density")
- `legend_shrink` (float): Fraction of original size to shrink the colorbar (default: 0.8)
- `**kwargs`: Additional arguments for plotting

#### Returns
- `matplotlib.axes.Axes` or `list of matplotlib.axes.Axes`: Axes object(s) with the plot(s)

#### Examples

```python
# Plot density for a single gene
pynebulosa.plot_density(adata, 'gene_0')

# Plot density for multiple genes
pynebulosa.plot_density(adata, ['gene_0', 'gene_1', 'gene_2'])

# Plot joint density for multiple genes
pynebulosa.plot_density(adata, ['gene_0', 'gene_1'], joint=True)

# Customize appearance
pynebulosa.plot_density(adata, 'gene_0', figsize=(10, 8), legend_title="Expression Level")

# Hide legend title
pynebulosa.plot_density(adata, 'gene_0', show_legend_title=False)
```

## Internal Functions

### `_plot_final_density`

Internal function that handles the actual plotting logic after data preprocessing.

### `plot_density_`

Lowest level plotting function that creates the actual matplotlib plots.

### `_extract_feature_data`

Extracts feature data from an AnnData object.

### `_search_dimensions`

Validates and selects dimensions from embeddings.

### `_validate_dimensions`

Validates that exactly two dimensions are provided.

## Utility Functions

### `calculate_density`

Calculates kernel density estimates for gene expression data.

## Parameter Details

### `show_legend_title`
Controls whether the legend title is displayed:
- Default: True (shows legend title)
- Hidden: False (hides legend title)

### `legend_title`
Controls the text displayed on the colorbar legend:
- Default: "Density"
- Custom text: Any string
- Used only when `show_legend_title=True`

### `legend_shrink`
Controls the size of the colorbar legend:
- Range: 0.0 to 1.0
- Default: 0.8
- Full size: 1.0
- Half size: 0.5

### `figsize`
Controls the figure size in inches:
- Format: Tuple of (width, height)
- Default: None (uses matplotlib defaults)
- Example: (10, 8) for a 10x8 inch figure

### `basis`
Specifies which embedding to use for visualization:
- Default: Last embedding in adata.obsm
- Examples: 'umap', 'pca', 'tsne'

### `dims`
Selects which dimensions to plot:
- Format: List of two integers
- Default: [1, 2] (first and second dimensions)
- Example: [1, 3] to plot first and third dimensions

## Best Practices

1. **Performance**: For large datasets, consider setting `raster=True` to improve rendering performance
2. **Visualization**: Use `joint=True` to visualize co-expression patterns between multiple genes
3. **Customization**: Adjust `adjust` parameter to control smoothing (values < 1 increase smoothing, values > 1 decrease smoothing)
4. **Color Palettes**: Use any matplotlib colormap name for the `pal` parameter
## See Also

- [AnnData Documentation](https://anndata.readthedocs.io/)
- [Matplotlib Documentation](https://matplotlib.org/)

- [Scanpy Documentation](https://scanpy.readthedocs.io/)
