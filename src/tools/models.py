import numpy as np
from numpy import linalg

x = 0
y = 1
z = 2

def base_model(mu, t, X):
    """
    Computes the next time derivative step to be used in an ode solver. It uses
    a basic gravitaional model that only takes the gravity of the central body
    into account.

    Parameters
    ----------
    mu : The gravitational parameter of the central body
    t  : The current time step
    X  : The current state vector

    Returns
    -------
    dX : The calculated values of the next time derivative step
    """
    r = np.array(X[:3]).reshape((3,1))
    v = np.array(X[3:6]).reshape((3,1))

    a = -mu/linalg.norm(r)**3*np.eye(3)@r

    dX = np.zeros(6)
    dX[:3] = v.reshape((1,3))
    dX[3:6] = a.reshape((1,3))
    return dX

def jgm3_model(mu, J2, J3, t, X):
    """
    Computes the next time derivative step to be used in an ode solver. It uses
    a model based on the joint gravity model 3, which includes the gravity of
    the central body and the second and third zonal harmonic.

    Parameters
    ----------
    mu : The gravitational parameter of the central body
    J2 : The second zonal harmonic
    J3 : The third zonal harmonic
    t  : The current time step
    X  : The current state vector

    Returns
    -------
    dX : The calculated values of the next time derivative step
    """
    r = np.array(X[:3]).reshape((3,1))
    v = np.array(X[3:6]).reshape((3,1))

    a = -mu/linalg.norm(r)**3*np.eye(3)@r
    a -= 3/2*mu/linalg.norm(r)**5*R**2*J2*(5*(r[z]/linalg.norm(r))**2*np.eye(3) - np.array([[1,0,0],[0,1,0],[0,0,3]]))@r
    a -= 1/2*mu/linalg.norm(r)**5*R**5*J3*(35*(r[z]/linalg.norm(r)**3 - 15*np.array([[1,0,0],[0,1,0],[0,0,2]])*r[z]/linalg.norm(r) + np.array([[0,0,0],[0,0,0],[0,0,3]])))@r

    dX = np.zeros(6)
    dX[:3] = v.reshape((1,3))
    dX[3:6] = a.reshape((1,3))
    return dX

def CR3BP(mu, t, X):
    """
    Computes the next time derivative step to be used in an ode solver. It uses
    a non-dimenial circular restricted 3-body problem model.

    Parameters
    ----------
    mu : The non-dimensional gravitational parameter for the system
    t  : The current time step
    X  : The current state vector

    Returns
    -------
    dX : The calculated values of the next time derivative step
    """
    r = np.array(X[:3]).reshape((3,1))
    v = np.array(X[:3]).reshape((3,1))

    r1 = np.sqrt((r[x] + mu)**2 + r[y]**2 + r[z]**2)
    r2 = np.sqrt((r[x] - 1 + mu)**2 + r[y]**2 + r[z]**2)

    dX = np.zeros(6)
    dX[:3] = v.reshape((1,3))
    dX[3] = r[x] + 2*v[y] - (1 - mu)*(r[x] + mu)/r1**3 - mu*(r[x] - 1 + mu)/r2**3
    dX[4] = r[y] - 2*v[x] - (1 - mu)*r[y]/r1**3 - mu*r[y]/r2**3
    dX[5] = -(1 - mu)*r[z]/r1**3 - mu*r[z]/r2**3
