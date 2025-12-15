"""
Plotting functions for PyNebulosa
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize
import anndata
from typing import Union, List, Optional, Dict, Any, Tuple
from .kde import calculate_density


def plot_density_(z: np.ndarray, 
                  feature: str, 
                  cell_embeddings: np.ndarray, 
                  dim_names: List[str], 
                  ax: Optional[plt.Axes] = None,
                  basis: Optional[str] = None,
                  shape: str = 'o', 
                  size: float = 1, 
                  pal: str = "viridis",
                  show_legend_title: bool = True,
                  legend_title: str = "Density",
                  legend_shrink: float = 0.8,
                  figsize: Optional[Tuple[float, float]] = None,
                  raster: bool = True,
                  **kwargs) -> plt.Axes:
    """
    Plot density estimates
    
    Parameters
    ----------
    z : array-like
        Vector with density values for each cell. Higher values indicate higher density regions.
    feature : str
        Name of the feature being plotted (e.g., gene name). Used as the plot title.
    cell_embeddings : array-like
        Matrix with cell embeddings of shape (n_cells, n_dimensions). 
        Typically UMAP, PCA, or t-SNE coordinates.
    dim_names : list
        Names of the dimensions from the cell embeddings. Used for axis labels.
    ax : matplotlib.axes.Axes, optional
        Pre-existing axes object to plot on. If None, a new figure and axes are created.
    basis : str, optional
        Name of the embedding basis (e.g., 'umap', 'pca'). 
        Used to create proper axis labels (e.g., UMAP1, UMAP2).
    shape : str
        Point shape for the scatter plot. Follows matplotlib marker conventions.
        Default: 'o' (circle)
    size : float
        Size of the points in the scatter plot. Multiplied by 10 for actual display size.
        Default: 1
    pal : str
        Color palette name for the density visualization. Any matplotlib colormap name.
        Default: "viridis"
    show_legend_title : bool
        Whether to show legend title on the colorbar. 
        When False, the colorbar legend title is completely hidden.
        Default: True
    legend_title : str
        String used as legend title on the colorbar. 
        Only displayed when show_legend_title=True.
        Default: "Density"
    legend_shrink : float
        Fraction of original size to shrink the colorbar. 
        Value between 0.0 (invisible) and 1.0 (full size).
        Default: 0.8
    figsize : tuple of float, optional
        Figure size as (width, height) in inches. 
        If None, uses matplotlib default (6.4, 4.8).
        Default: None
    raster : bool
        Whether to rasterize the plot for better performance with large datasets.
        Recommended for datasets with >1000 cells.
        Default: True
    **kwargs : dict
        Additional arguments passed to matplotlib scatter plot function.
        
    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot
    """
    # Create DataFrame for plotting
    plot_data = pd.DataFrame({
        dim_names[0]: cell_embeddings[:, 0],
        dim_names[1]: cell_embeddings[:, 1],
        'feature': z
    })
    
    # Create axes if not provided
    if ax is None:
        # Use custom figsize if provided, otherwise use default
        if figsize is not None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig, ax = plt.subplots(figsize=(6.4, 4.8))  # Matplotlib default figure size
    
    # Create scatter plot
    scatter = ax.scatter(plot_data[dim_names[0]], 
                        plot_data[dim_names[1]], 
                        c=plot_data['feature'],
                        s=size*10,
                        cmap=pal,
                        edgecolors='none',
                        rasterized=raster,
                        **kwargs)
    
    # Labels and title (following Scanpy style)
    # For embedding plots, use the basis name to create proper labels (e.g., UMAP1, UMAP2)
    if basis is not None:
        basis_clean = basis.replace('X_', '').upper()
        xlabel = f"{basis_clean}1"
        ylabel = f"{basis_clean}2"
    else:
        # Fallback to generic dimension names
        xlabel = dim_names[0]
        ylabel = dim_names[1]
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(feature)
    
    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=legend_shrink)
    if show_legend_title and legend_title != "":
        cbar.set_label(legend_title)
    
    # Styling to match Scanpy defaults
    # Keep all spines visible (Scanpy default)
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)

    # Turn off grid
    ax.grid(False)
    # Set spine linewidth to match Scanpy defaults
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
    
    return ax


def _validate_dimensions(dims):
    """
    Validate that exactly two dimensions are provided
    
    Parameters
    ----------
    dims : list of int
        List of dimension indices to validate. Should contain exactly two integers.
        
    Raises
    ------
    ValueError
        If dims does not contain exactly two dimensions.
    """
    if len(dims) != 2:
        raise ValueError("Only two dimensions can be plotted")


def _search_dimensions(dims, cell_embeddings, reduction):
    """
    Validate and select dimensions from embeddings
    
    Parameters
    ----------
    dims : list of int
        List of dimension indices to select (1-based indexing).
    cell_embeddings : np.ndarray
        Matrix with cell embeddings of shape (n_cells, n_dimensions).
    reduction : str
        Name of the dimensionality reduction technique (e.g., 'umap', 'pca').
        
    Returns
    -------
    np.ndarray
        Selected dimensions from cell_embeddings as a 2D array of shape (n_cells, 2).
        
    Raises
    ------
    ValueError
        If any dimension in dims is not present in cell_embeddings.
    """
    n_dims = cell_embeddings.shape[1] if cell_embeddings.ndim > 1 else 1
    available_dims = list(range(1, n_dims + 1)) if n_dims > 0 else []
    
    invalid_dims = [d for d in dims if d not in available_dims]
    if invalid_dims:
        raise ValueError(f"Dimension(s) {invalid_dims} not present in {reduction}")
    
    # Convert to 0-based indexing
    idx_dims = [d - 1 for d in dims]
    
    if cell_embeddings.ndim == 1:
        # If only one dimension, return as column vector
        return cell_embeddings.reshape(-1, 1)
    else:
        return cell_embeddings[:, idx_dims]


def _extract_feature_data(adata, features, layer=None):
    """
    Extract feature data from AnnData object
    
    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object containing expression data.
    features : list of str
        List of feature names to extract.
    layer : str, optional
        Layer to use for expression data. If None, uses .X (main data matrix).
        
    Returns
    -------
    pd.DataFrame
        DataFrame with feature data of shape (n_cells, n_features).
        Index contains cell names, columns contain feature names.
        
    Raises
    ------
    ValueError
        If any feature in features is not present in adata.var_names.
    """
    # Check if features exist in var_names
    missing_features = [f for f in features if f not in adata.var_names]
    if missing_features:
        raise ValueError(f"'{', '.join(missing_features)}' feature(s) not present in var_names")
    
    # Extract data
    if layer is None:
        # Use X (default)
        data = adata[:, features].X
    else:
        if layer not in adata.layers:
            raise ValueError(f"Layer '{layer}' not found in AnnData object")
        data = adata[:, features].layers[layer]
    
    # Convert to dense array if sparse
    if hasattr(data, 'toarray'):
        data = data.toarray()
    
    # Ensure it's a 2D array
    if data.ndim == 1:
        data = data.reshape(-1, 1)
    
    # Create DataFrame with cell names as index and feature names as columns
    df = pd.DataFrame(data, index=adata.obs_names, columns=features)
    return df


def _plot_final_density(vars_df: pd.DataFrame, 
                       cell_embeddings: np.ndarray, 
                       features: List[str], 
                       ax: Optional[plt.Axes] = None,
                       basis: Optional[str] = None,
                       joint: bool = False,
                       only_joint_plot: bool = False,
                       method: str = 'wkde',
                       adjust: float = 1,
                       shape: str = 'o',
                       size: float = 1,
                       pal: str = "viridis",
                       raster: bool = True,
                       figsize: Optional[Tuple[float, float]] = None,
                       show_legend_title: bool = True,
                       legend_title: str = "Density",
                       legend_shrink: float = 0.8,
                       **kwargs) -> Union[plt.Axes, List[plt.Axes]]:
    """
    Plot final density
    
    Parameters
    ----------
    vars_df : pd.DataFrame
        DataFrame with feature data of shape (n_cells, n_features).
        Index contains cell names, columns contain feature names.
    cell_embeddings : np.ndarray
        Matrix with cell embeddings of shape (n_cells, n_dimensions).
        Typically UMAP, PCA, or t-SNE coordinates.
    features : list of str
        List of features to plot (e.g., gene names).
    ax : matplotlib.axes.Axes, optional
        Pre-existing axes object to plot on. If None, new figure(s) and axes are created.
    basis : str, optional
        Name of the embedding basis (e.g., 'umap', 'pca').
        Used to create proper axis labels (e.g., UMAP1, UMAP2).
    joint : bool
        Return joint density plot? When True, displays both individual gene densities
        and joint density in subplots. For single features, this parameter has no effect.
        Default: False
    only_joint_plot : bool
        When True and joint=True, only plot the joint density map (ignore individual gene plots).
        This suppresses individual gene density plots and shows only the combined joint density.
        Default: False
    method : str
        Kernel density estimation method. Currently only 'wkde' (weighted kernel density estimation) is supported.
        Default: 'wkde'
    adjust : float
        Numeric value to adjust the bandwidth of the kernel density estimation.
        Values < 1.0 increase smoothing, values > 1.0 decrease smoothing.
        Default: 1
    shape : str
        Point shape for the scatter plot. Follows matplotlib marker conventions.
        Default: 'o' (circle)
    size : float
        Size of the points in the scatter plot. Multiplied by 10 for actual display size.
        Default: 1
    pal : str
        Color palette name for the density visualization. Any matplotlib colormap name.
        Default: "viridis"
    raster : bool
        Whether to rasterize the plot for better performance with large datasets.
        Recommended for datasets with >1000 cells.
        Default: True
    figsize : tuple of float, optional
        Figure size as (width, height) in inches.
        For multiple subplots, this size is multiplied by the grid dimensions.
        If None, uses matplotlib default (6.4, 4.8).
        Default: None
    show_legend_title : bool
        Whether to show legend title on the colorbar.
        When False, the colorbar legend title is completely hidden.
        Default: True
    legend_title : str
        String used as legend title on the colorbar.
        Only displayed when show_legend_title=True.
        For joint plots, individual gene plots use this title, but joint density plots use "Joint density".
        Default: "Density"
    legend_shrink : float
        Fraction of original size to shrink the colorbar.
        Value between 0.0 (invisible) and 1.0 (full size).
        Default: 0.8
    **kwargs : dict
        Additional arguments passed to matplotlib scatter plot function.
        
    Returns
    -------
    matplotlib.axes.Axes or list of matplotlib.axes.Axes
        Axes object(s) with the plot(s). For single plots, returns a single Axes object.
        For multiple plots, returns a numpy array of Axes objects.
    """
    dim_names = [f"dim_{i}" for i in range(cell_embeddings.shape[1])] if cell_embeddings.shape[1] <= 2 else [f"dim_0", f"dim_1"]
    
    # When joint=True, always create combined plots
    # (Previously this would set combine=True, but we're removing the combine parameter)
    
    if vars_df.shape[1] > 1:
        # Multiple features
        res_list = []
        axes_list = []
        
        # When only_joint_plot=True and joint=True, plot only the joint density
        if only_joint_plot and joint:
            # Calculate joint density
            joint_z = np.ones(len(vars_df))
            for feature in vars_df.columns:
                z = calculate_density(vars_df[feature].values, cell_embeddings, method, adjust)
                joint_z *= z
            
            # Create single plot
            joint_label = " + ".join(features)
            # For only joint plot, use "Joint density" as legend title when show_legend_title=True
            joint_legend_title = "Joint density" if show_legend_title else ""
            return plot_density_(joint_z, joint_label, cell_embeddings, dim_names, ax, basis, shape, size, 
                               pal, show_legend_title, joint_legend_title, legend_shrink, figsize, raster, **kwargs)
        elif joint:
            # When joint=True, plot individual gene densities and joint density
            # Create subplots
            n_features = len(vars_df.columns) + 1  # +1 for joint plot
            n_cols = min(3, n_features)
            n_rows = (n_features + n_cols - 1) // n_cols
            
            # Calculate figure size
            if figsize is not None:
                fig_width = figsize[0] * n_cols
                fig_height = figsize[1] * n_rows
            else:
                fig_width = 6.4 * n_cols
                fig_height = 4.8 * n_rows
                
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))
            if n_features > 1:
                axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]
            else:
                axes = [axes]
                
            # Plot individual gene densities
            for i, feature in enumerate(vars_df.columns):
                z = calculate_density(vars_df[feature].values, cell_embeddings, method, adjust)
                plot_density_(z, feature, cell_embeddings, dim_names, axes[i], basis, shape, size, 
                             pal, show_legend_title, legend_title, legend_shrink, figsize, raster, **kwargs)
            
            # Calculate and plot joint density
            joint_z = np.ones(len(vars_df))
            for feature in vars_df.columns:
                z = calculate_density(vars_df[feature].values, cell_embeddings, method, adjust)
                joint_z *= z
            
            joint_label = " + ".join(features)
            # For joint plot, use "Joint density" as legend title when show_legend_title=True
            joint_legend_title = "Joint density" if show_legend_title else ""
            plot_density_(joint_z, joint_label, cell_embeddings, dim_names, axes[len(vars_df.columns)], basis, shape, size, 
                         pal, show_legend_title, joint_legend_title, legend_shrink, figsize, raster, **kwargs)
            
            # Hide unused subplots
            for i in range(len(vars_df.columns) + 1, len(axes)):
                axes[i].set_visible(False)
                
            # Apply tight_layout to prevent overlapping labels
            plt.tight_layout()
                
            return axes[0] if len(axes) == 1 else axes
        else:
            # When joint=False, plot individual gene densities in subplots
            # (This replaces the previous "combine=True" behavior)
            # For combined plot, we'll handle this differently
            for feature in vars_df.columns:
                z = calculate_density(vars_df[feature].values, cell_embeddings, method, adjust)
                res_list.append(z)
            
            # Create subplots
            n_features = len(vars_df.columns)
            n_cols = min(3, n_features)
            n_rows = (n_features + n_cols - 1) // n_cols
            
            # Calculate figure size
            if figsize is not None:
                fig_width = figsize[0] * n_cols
                fig_height = figsize[1] * n_rows
            else:
                fig_width = 6.4 * n_cols
                fig_height = 4.8 * n_rows
                
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))
            if n_features > 1:
                axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]
            else:
                axes = [axes]
                
            # Plot individual features
            for i, feature in enumerate(vars_df.columns):
                z = calculate_density(vars_df[feature].values, cell_embeddings, method, adjust)
                plot_density_(z, feature, cell_embeddings, dim_names, axes[i], basis, shape, size, 
                             pal, show_legend_title, legend_title, legend_shrink, figsize, raster, **kwargs)
            
            # Hide unused subplots
            for i in range(len(vars_df.columns), len(axes)):
                axes[i].set_visible(False)
                
            # Apply tight_layout to prevent overlapping labels
            plt.tight_layout()
                
            return axes[0] if len(axes) == 1 else axes
    else:
        # Single feature
        z = calculate_density(vars_df.iloc[:, 0].values, cell_embeddings, method, adjust)
        return plot_density_(z, features[0], cell_embeddings, dim_names, ax, basis, shape, size, 
                           pal, show_legend_title, legend_title, legend_shrink, figsize, raster, **kwargs)


def plot_density(adata: anndata.AnnData,
                 features: Union[str, List[str]],
                 ax: Optional[plt.Axes] = None,
                 basis: Optional[str] = None,
                 layer: Optional[str] = None,
                 dims: List[int] = [1, 2],
                 joint: bool = False,
                 only_joint_plot: bool = False,
                 method: str = 'wkde',
                 adjust: float = 1,
                 shape: Union[str, int] = 'o',
                 size: float = 1,
                 pal: str = "viridis",
                 raster: bool = True,
                 figsize: Optional[Tuple[float, float]] = None,
                 show_legend_title: bool = True,
                 legend_title: str = "Density",
                 legend_shrink: float = 0.8,
                 **kwargs) -> Union[plt.Axes, List[plt.Axes]]:
    """
    Plot gene-weighted 2D kernel density for AnnData objects
    
    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object with cell embeddings and gene expression data.
        Must contain .X (or specified layer) with expression data and .obsm with embeddings.
    features : str or list of str
        Features (e.g. genes) to visualize. Can be a single feature name or list of feature names.
    ax : matplotlib.axes.Axes, optional
        Pre-existing axes object to plot on. If None, new figure(s) and axes are created.
    basis : str, optional
        Name of the embedding to visualize (e.g., 'umap', 'pca').
        If not provided, uses the last embedding in .obsm.
        Default: None
    layer : str, optional
        Layer to use for expression data. If None, uses .X (main data matrix).
        Must be a valid layer name in adata.layers.
        Default: None
    dims : list of int
        Vector of length 2 specifying the dimensions to be plotted.
        Uses 1-based indexing (e.g., [1, 2] for first and second dimensions).
        Default: [1, 2]
    joint : bool
        Return joint density plot? When True, displays both individual gene densities
        and joint density in subplots. For single features, this parameter has no effect.
        Default: False
    only_joint_plot : bool
        When True and joint=True, only plot the joint density map (ignore individual gene plots).
        This suppresses individual gene density plots and shows only the combined joint density.
        Default: False
    method : str
        Kernel density estimation method. Currently only 'wkde' (weighted kernel density estimation) is supported.
        Default: 'wkde'
    adjust : float
        Numeric value to adjust the bandwidth of the kernel density estimation.
        Values < 1.0 increase smoothing, values > 1.0 decrease smoothing.
        Default: 1
    shape : str or int
        Point shape for the scatter plot. Follows matplotlib marker conventions.
        Default: 'o' (circle)
    size : float
        Size of the points in the scatter plot. Multiplied by 10 for actual display size.
        Default: 1
    pal : str
        Color palette name for the density visualization. Any matplotlib colormap name.(e.g., 'viridis', 'plasma', 'inferno', 'magma', 'Blues', 'Reds','RdBu_r' etc.)
        Default: "viridis"
    raster : bool
        Whether to rasterize the plot for better performance with large datasets.
        Recommended for datasets with >1000 cells.
        Default: True
    figsize : tuple of float, optional
        Figure size as (width, height) in inches.
        For multiple subplots, this size is multiplied by the grid dimensions.
        If None, uses matplotlib default (6.4, 4.8).
        Default: None
    show_legend_title : bool
        Whether to show legend title on the colorbar.
        When False, the colorbar legend title is completely hidden.
        Default: True
    legend_title : str
        String used as legend title on the colorbar.
        Only displayed when show_legend_title=True.
        For joint plots, individual gene plots use this title, but joint density plots use "Joint density".
        Default: "Density"
    legend_shrink : float
        Fraction of original size to shrink the colorbar.
        Value between 0.0 (invisible) and 1.0 (full size).
        Default: 0.8
    **kwargs : dict
        Additional arguments passed to matplotlib scatter plot function.
        
    Returns
    -------
    matplotlib.axes.Axes or list of matplotlib.axes.Axes
        Axes object(s) with the plot(s). For single plots, returns a single Axes object.
        For multiple plots, returns a numpy array of Axes objects.
    """
    # Validate dimensions
    _validate_dimensions(dims)
    
    # Handle single feature
    if isinstance(features, str):
        features = [features]
    
    # Set up default basis
    if basis is None:
        # Use the last embedding in obsm
        if len(adata.obsm) == 0:
            raise ValueError("No embeddings found in AnnData object")
        basis = list(adata.obsm.keys())[-1]
    
    # Validate basis exists
    basis_key = basis if basis.startswith('X_') else f'X_{basis}'
    if basis_key not in adata.obsm:
        # Try without prefix
        if basis not in adata.obsm:
            raise ValueError(f"No embedding named '{basis}' found in AnnData object")
        basis_key = basis
    
    # Get cell embeddings
    cell_embeddings = adata.obsm[basis_key]
    
    # Search for dimensions
    cell_embeddings = _search_dimensions(dims, cell_embeddings, basis)
    
    # Extract feature data
    vars_df = _extract_feature_data(adata, features, layer)
    
    # Plot final density
    return _plot_final_density(vars_df, cell_embeddings, features, ax, basis, joint, only_joint_plot,
                              method, adjust, shape, size, pal,
                              raster, figsize, show_legend_title, legend_title, legend_shrink, **kwargs)

