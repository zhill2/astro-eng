import numpy as np

# ROTATIONS
def xrot(omega):
    """
    Computes the rotation matrix about the x-axis

    Parameters
    ----------
    omega : The rotation angle in radians

    Returns
    -------
    The rotation matrix
    """
    return np.array([[1, 0, 0], [0, np.cos(omega), np.sin(omega)],[0, -np.sin(omega), np.cos(omega)]]);

def yrot(omega):
    """
    Computes the rotation matrix about the y-axis

    Parameters
    ----------
    omega : The rotation angle in radians

    Returns
    -------
    The rotation matrix
    """
    return np.array([[np.cos(omega), 0, -np.sin(omega)], [0, 1, 0], [np.sin(omega), 0, np.cos(omega)]]);

def zrot(omega):
    """
    Computes the rotation matrix about the z-axis

    Parameters
    ----------
    omega : The rotation angle in radians

    Returns
    -------
    The rotation matrix
    """
    return np.array([[np.cos(omega), np.sin(omega), 0],[-np.sin(omega), np.cos(omega), 0], [0, 0, 1]]);
