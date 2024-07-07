import numpy as np
from . import transforms

tle_format = {
    'sat_num': [1,3,7],
    'class': [1,8,8],
    'ID_year': [1,10,11],
    'ID_launch': [1,12,14],
    'ID_piece': [1,15,17],
    'epoch_year': [1,19,20],
    'epoch_day': [1,21,32],
    'dMdt': [1,34,43],
    'd2Mdt2': [1,45,52],
    'bstar': [1,54,61],
    'eph_type': [1,63,63],
    'elem_num': [1,65,68],
    'i': [2,9,16],
    'RAAN': [2,18,25],
    'e': [2,27,33],
    'omega': [2,35,42],
    'M': [2,44,51],
    'n': [2,53,63],
    'rev': [2,64,68]
}

def parsetle(tle):
    """
    Parses a given two-line element.

    Parameters
    ----------
    tle : A list containing the TLE

    Returns
    -------
    parsed : A dictionary conatining the complete parsed TLE
    """"
    parsed = {}
    parsed['name'] = tle[0]
    for key, val in tle_format.items():
        line = val[0]
        start = val[1] - 1
        end = val[2]
        if (end - start) > 0:
            parsed[key] = tle[line][start:end]
        else:
            parsed[key] = tle[line][start]
    return parsed

def tle2coes(tle, mu):
    """
    Pulls the classical orbital elements out of the given TLE.

    Parameters
    ----------
    tle : A list containing the TLE
    mu  : The gravitaional parameter of the central body

    Returns
    -------
    a     : The semimajor axis of the orbit
    i     : The inclination of the orbit
    RAAN  : The right ascension of the ascending node
    e     : The eccentricity of the orbit
    omega : The argument of periapsis of the orbit
    M     : The mean motion of the orbit
    """
    parsed = parsetle(tle)
    
    i = float(parsed['i'])
    RAAN = float(parsed['RAAN'])
    e = float("0." + parsed['e'])
    omega = float(parsed['omega'])
    M = float(parsed['M'])
    n = float(parsed['n'])
    a = (mu/(n*2*np.pi/(3600*24))**2)**(1/3)

    return [a, i, RAAN, e, omega, M]

def tle2state(tle, mu):
    """
    Pulls the classical orbital elements out of the given TLE.

    Parameters
    ----------
    tle : A list containing the TLE
    mu  : The gravitaional parameter of the central body

    Returns
    -------
    Returns the state of the spacecraft at the time of the TLE
    """
    coes = tle2coes(tle, mu)
    return transforms.coes2state(coes,mu)
