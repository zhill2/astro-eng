import numpy as np
from . import matrix

def car2sph(x,y,z):
    r = np.sqrt(x**2 + y**2 + z**2)
    if (x == 0) and (y == 0):
        long = 0.0
        if z < 0:
            lat = -90.0
        else:
            lat = 90.0
    else:
        long = np.arctan2(y,x)*180.0/np.pi
        lat = np.arctan2(sqrt(x**2 + y**2),z)*180.0/np.pi
    return [r, lat, long]

def sph2car(r, lat, long):
    colat = np.pi/2.0 - lat*np.pi/180.0
    x = r*np.sin(colat)*np.cos(long*np.pi/180.0)
    y = r*np.sin(colat)*np.sin(long*np.pi/180.0)
    z = r*np.cos(colat)
    return [x, y, z]

def eci2ecef(x, y, z, t):
    pass
    

def ecef2eci(r, lat, long, t):
    pass

def coes2state(coes, mu):
    a = coes[0]
    i = coes[1]
    RAAN = coes[2]
    e = coes[3]
    omega = coes[4]
    M = coes[5]
    h = np.sqrt(mu*a*(1 - e)**2)

    diff = 1
    E = M
    while abs(diff) > 1e-6:
        E_new = M + e*np.sin(E)
        diff = E - E_new
        E = E_new
        
    nu = 2*np.arctan(np.tan(E/2)/np.sqrt((1 - e)/(1 + e)))
    r = h**2/mu/(1 + e*np.cos(nu))*np.array([np.cos(nu), np.sin(nu), 0])
    v =mu/h*np.array([-np.sin(nu), e + np.cos(nu), 0])

    rot1 = matrix.zrot(omega)
    rot2 = matrix.xrot(i)
    rot3 = matrix.zrot(RAAN)

    r = r@rot1@rot2@rot3
    v = v@rot1@rot2@rot3

    return np.concatenate((r,v))

def state2coes(state, mu):
    pass
