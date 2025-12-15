"""
Kernel density estimation functions for PyNebulosa
"""

import numpy as np
from scipy.stats import norm
from sklearn.neighbors import KernelDensity
import warnings


def wkde2d(x, y, w, adjust=1, n=100, lims=None):
    """
    Weighted 2D kernel density estimation
    
    Parameters
    ----------
    x : array-like
        Dimension 1
    y : array-like
        Dimension 2
    w : array-like
        Weight variable
    adjust : float
        Bandwidth adjustment
    n : int
        Number of grid points in each direction
    lims : array-like
        The limits of the rectangle covered by the grid as [xl, xu, yl, yu]
        
    Returns
    -------
    dict
        A dictionary with keys 'x', 'y', 'z' containing the grid points and density estimates
    """
    # Convert to numpy arrays
    x = np.asarray(x)
    y = np.asarray(y)
    w = np.asarray(w)
    
    # Validate values and dimensions
    nx = len(x)
    if not (len(y) == nx and len(w) == nx):
        raise ValueError("data vectors must be the same length")
        
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        raise ValueError("missing or infinite values in the data are not allowed")
        
    if lims is None:
        lims = [np.min(x), np.max(x), np.min(y), np.max(y)]
    else:
        lims = np.asarray(lims)
        if not np.all(np.isfinite(lims)):
            raise ValueError("only finite values are allowed in 'lims'")
            
    # Bandwidth selection (using Silverman's rule of thumb as approximation)
    # For more accurate bandwidth selection, we could use a more sophisticated method
    def hpi(data):
        """Simple bandwidth selector (Silverman's rule of thumb)"""
        n = len(data)
        sigma = np.std(data, ddof=1)
        # Interquartile range
        q75, q25 = np.percentile(data, [75, 25])
        iqr = q75 - q25
        # Use smaller of std and IQR/1.349 as measure of spread
        sig = min(sigma, iqr/1.349) if iqr > 0 else sigma
        # Silverman's rule of thumb
        return sig * (4.0 / (3 * n)) ** (1/5)
    
    h = np.array([hpi(x), hpi(y)])
    h = h * adjust
    
    # Get grid
    gx = np.linspace(lims[0], lims[1], n)
    gy = np.linspace(lims[2], lims[3], n)
    
    # Compute density
    # Weight calculation
    ax = np.subtract.outer(gx, x) / h[0]
    ay = np.subtract.outer(gy, y) / h[1]
    
    # Create weight matrix
    w_matrix = np.tile(w, (n, 1))
    
    # Compute density using Gaussian kernel
    z = np.dot(norm.pdf(ax) * w_matrix, (norm.pdf(ay) * w_matrix).T) / (np.sum(w) * h[0] * h[1])
    
    return {'x': gx, 'y': gy, 'z': z}


def get_dens(data, dens_dict, method='wkde'):
    """
    Map density estimates to individual observations
    
    Parameters
    ----------
    data : array-like
        Data points to get density for
    dens_dict : dict
        Dictionary with density estimates
    method : str
        Method used for density estimation
        
    Returns
    -------
    array
        Density values for each observation
    """
    data = np.asarray(data)
    
    if method == 'wkde':
        # Find intervals for each data point
        x_indices = np.searchsorted(dens_dict['x'], data[:, 0], side='right') - 1
        y_indices = np.searchsorted(dens_dict['y'], data[:, 1], side='right') - 1
        
        # Clip indices to valid range
        x_indices = np.clip(x_indices, 0, len(dens_dict['x']) - 2)
        y_indices = np.clip(y_indices, 0, len(dens_dict['y']) - 2)
        
        # Get density values
        z_values = dens_dict['z'][x_indices, y_indices]
        
    else:
        raise ValueError(f"Method {method} not supported")
        
    return z_values


def calculate_density(w, x, method='wkde', adjust=1, map=True):
    """
    Estimate weighted kernel density
    
    Parameters
    ----------
    w : array-like
        Vector with weights for each observation (e.g., gene expression values)
    x : array-like
        Matrix with dimensions where to calculate the density from. Only
        the first two dimensions will be used
    method : str
        Kernel density estimation method ('wkde')
    adjust : float
        Numeric value to adjust to bandwidth
    map : bool
        Whether to map densities to individual observations
        
    Returns
    -------
    array or dict
        If map is True, a vector with corresponding densities for each observation
        is returned. Otherwise, a dictionary with the density estimates is returned.
    """
    w = np.asarray(w)
    x = np.asarray(x)
    
    # Ensure x is 2D
    if x.ndim == 1:
        x = x.reshape(-1, 1)
        
    # Use only first two dimensions
    if x.shape[1] < 2:
        raise ValueError("Need at least 2 dimensions in x")
        
    x_coords = x[:, 0]
    y_coords = x[:, 1]
    
    if method == 'wkde':
        dens = wkde2d(
            x=x_coords,
            y=y_coords,
            w=w / np.sum(w) * len(w),
            adjust=adjust
        )
    else:
        raise ValueError(f"Method {method} not supported. Only 'wkde' is currently supported.")
        
    if map:
        # Create data array for mapping
        data_points = np.column_stack([x_coords, y_coords])
        return get_dens(data_points, dens, method)
    else:
        return dens