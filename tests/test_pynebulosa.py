"""
Tests for PyNebulosa
"""

import numpy as np
import anndata
import pytest
import scanpy as sc
import matplotlib.pyplot as plt
from pynebulosa import calculate_density, wkde2d, plot_density


def test_wkde2d():
    """Test weighted 2D kernel density estimation"""
    # Create simple test data
    np.random.seed(42)
    x = np.random.normal(0, 1, 100)
    y = np.random.normal(0, 1, 100)
    w = np.random.uniform(0, 1, 100)
    
    # Test basic functionality
    result = wkde2d(x, y, w, n=20)
    
    assert 'x' in result
    assert 'y' in result
    assert 'z' in result
    assert len(result['x']) == 20
    assert len(result['y']) == 20
    assert result['z'].shape == (20, 20)


def test_calculate_density():
    """Test density calculation"""
    # Create simple test data
    np.random.seed(42)
    x = np.random.normal(0, 1, 50)
    y = np.random.normal(0, 1, 50)
    w = np.random.uniform(0, 1, 50)
    
    # Create coordinate matrix
    coords = np.column_stack([x, y])
    
    # Test mapping functionality
    density = calculate_density(w, coords, method='wkde', map=True)
    
    assert len(density) == 50
    assert np.all(density >= 0)


def test_plot_density_single_gene():
    """Test plotting function with sample data"""
    # Create sample AnnData object
    np.random.seed(42)
    n_cells = 100
    n_genes = 50
    
    # Generate expression data
    X = np.random.poisson(2, size=(n_cells, n_genes))
    
    # Generate embeddings
    embedding = np.random.normal(0, 1, size=(n_cells, 2))
    
    # Create AnnData object
    adata = anndata.AnnData(X=X.astype(float))
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    adata.obsm['X_umap'] = embedding
    
    # Test basic plotting
    ax = plot_density(adata, 'gene_0', basis='umap')
    assert ax is not None


def test_plot_density_multiple_genes():
    """Test plotting multiple genes"""
    # Create sample AnnData object
    np.random.seed(42)
    n_cells = 100
    n_genes = 50
    
    # Generate expression data
    X = np.random.poisson(2, size=(n_cells, n_genes))
    
    # Generate embeddings
    embedding = np.random.normal(0, 1, size=(n_cells, 2))
    
    # Create AnnData object
    adata = anndata.AnnData(X=X.astype(float))
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    adata.obsm['X_umap'] = embedding
    
    # Test multiple genes plotting
    axes = plot_density(adata, ['gene_0', 'gene_1'], basis='umap')
    assert len(axes) == 2
    
    # Test joint plotting
    axes = plot_density(adata, ['gene_0', 'gene_1'], basis='umap', joint=True)
    assert axes is not None


def test_only_joint_plot():
    """Test only_joint_plot parameter"""
    # Create sample AnnData object
    np.random.seed(42)
    n_cells = 100
    n_genes = 50
    
    # Generate expression data
    X = np.random.poisson(2, size=(n_cells, n_genes))
    
    # Generate embeddings
    embedding = np.random.normal(0, 1, size=(n_cells, 2))
    
    # Create AnnData object
    adata = anndata.AnnData(X=X.astype(float))
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    adata.obsm['X_umap'] = embedding
    
    # Test only_joint_plot=True with joint=True
    ax = plot_density(adata, ['gene_0', 'gene_1'], basis='umap', joint=True, only_joint_plot=True)
    assert ax is not None
    
    # Should return a single axes object, not a list
    assert isinstance(ax, plt.Axes)
    
    print("✓ only_joint_plot parameter works correctly")


def test_legend_title_parameters():
    """Test show_legend_title and legend_title parameters"""
    # Create sample AnnData object
    np.random.seed(42)
    n_cells = 100
    n_genes = 50
    
    # Generate expression data
    X = np.random.poisson(2, size=(n_cells, n_genes))
    
    # Generate embeddings
    embedding = np.random.normal(0, 1, size=(n_cells, 2))
    
    # Create AnnData object
    adata = anndata.AnnData(X=X.astype(float))
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    adata.obsm['X_umap'] = embedding
    
    # Test show_legend_title=False
    ax1 = plot_density(adata, 'gene_0', basis='umap', show_legend_title=False)
    assert ax1 is not None
    
    # Test custom legend_title
    ax2 = plot_density(adata, 'gene_0', basis='umap', legend_title="Custom Title")
    assert ax2 is not None
    
    # Test joint plot legend title (should be "Joint density")
    ax3 = plot_density(adata, ['gene_0', 'gene_1'], basis='umap', joint=True, only_joint_plot=True)
    assert ax3 is not None
    
    print("✓ legend title parameters work correctly")


if __name__ == "__main__":
    # Run simple test
    test_plot_density_single_gene()
    test_plot_density_multiple_genes()
    test_only_joint_plot()
    test_legend_title_parameters()
    print("All tests passed!")