import numpy as np
import matplotlib.pyplot as plt

def create_surf(R):
    """
    Creates a spherical surface that can be used to represent a planet

    Parameter
    ---------
    R : Radius of the sphere

    Returns
    -------
    A matrix of x,y,z coordinates for a sphere of radius R
    """
    phi, theta = np.mgrid[0.0:np.pi:20j, 0.0:2.0*np.pi:20j]

    x = R*np.sin(phi)*np.cos(theta)
    y = R*np.sin(phi)*np.sin(theta)
    z = R*np.cos(phi)

    return [x,y,z]

def plot_orbit(states, surf=None):
    """
    Creates a plot of the given states and surface

    Parameter
    ---------
    states : A matrix of the states of the spacecraft in cartesian coordinates
    surf   : A matrix of x,y,z coordinates for a sphere

    Returns
    -------
    Nothing
    """
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot(states[0], states[1], states[2], color='b', zorder=10)
    if surf != None:
        ax.plot_wireframe(surf[0], surf[1], surf[2], color='g', zorder=0)
    plt.axis('equal')
    plt.title('Orbit')
    plt.show()

def plot_ground(states):
    """
    Creates a plot of the given lat and long coordinates of the spacecraft

    Parameter
    ---------
    states : A matrix of the states of the spacecraft in spherical coordinates

    Returns
    -------
    Nothing
    """
    plt.figure()
    ax = plt.axes()
    ax.plot(states[1], states[2])
    plt.axis('equal')
    plt.title('Lat-Long')
    plt.show()
