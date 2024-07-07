import context
import numpy as np
from tools import propogator

def cholesky(M,N):
    """
    Solves the problem Mx=N using the Cholesky decomposition method.

    Parameters
    ----------
    M : An nxn matrix
    N : An nx1 vector

    Returns
    -------
    x : The solution to the given problem
    P : The covariance matrix for the given problem

    Raises
    ------
    Exception : if the parameters are not of the right shape

    Example
    -------
    >>> import numpy as np
    >>> M = [[6.01, -3], [-3, 6.01]]
    >>> N = [3.12, 2.82]
    >>> [x, P] = cholesky(M,N)
    >>> print(x)
    """
    M = np.array(M)
    N = np.array(N)
    print(max(N.shape))

    if M.shape[0] != M.shape[1]:
        raise Exception("The M matrix must be square")
    if M.shape[0] != max(N.shape):
        raise Exception("The M matrix size must agree with the N vector")
    
    n = M.shape[0]

    R = np.zeros(M.shape)
    for i in range(n):
        calc = 0
        for k in range(i-1):
            calc += R[k][i]**2

        R[i][i] = np.sqrt(M[i][i] - calc)

        for j in range(i-1):
            calc = 0
            for k in range(i-1):
                calc += R[k][i]*R[k][j]

            R[i][j] = (M[i][j] - calc)/R[i][i]

    Z = np.zeros((n,1))
    for i in range(n):
        calc = 0
        for j in range(i-1):
            calc += R[j][i]*Z[j]
        Z[i] = (N[i] - calc)/R[i][i]

    xhat = np.zeros((n,1))
    for i in reversed(range(n)):
        calc = 0
        for j in range(i+1,n):
            calc += R[i][j]*xhat[j]
        xhat[i] = (Z[i] - calc)/R[i][i]

    S = np.zeros(M.shape)
    for i in range(n):
        S[i][i] = 1/R[i][i]
    for i in range(n):
        for j in range(n):
            calc = 0
            for k in range(i,j-1):
                calc += R[k][j]*S[i][k]
            S[i][j] = -S[i][j]*calc

    P = S*np.transpose(S)
    return [xhat, P]

def std_dev(P):
    sigma = np.zeros(len(P))
    for i in range(len(P)):
        sigma[i] = np.sqrt(P[i][i])
    return sigma

def ode(t,X):
    R = 6378136.3
    S = 3
    m = 970
    rho0 = 3.614e-13
    r0 = 700000.0
    H = 88667
    thetadot = 7.29211158543e-5

    x = X[i][0]
    y = X[i][1]
    z = X[i][2]
    xdot = X[i][3]
    ydot = X[i][4]
    zdot = X[i][5]
    mu = X[i][6]
    J2 = X[i][7]
    CD = X[i][8]
    phi = reshape(X[i][18:],18,18)

    r = np.sqrt(x**2 + y**2 + z**2)
    rho = rho0*np.exp(-(r - r0)/H)
    v_a = linalg.norm([xdot,ydot,zdot] - np.cross([0,0,thetadot],[x,y,z]))
    A = Amatrix(X)

    dX = np.zeros(18 + 18*18,1)
    dX[0] = xdot
    dX[1] = ydot
    dX[2] = zdot
    dX[3] = -mu*x*(1/r**3 + 3/2*J2*R**2/r**5*(1 - 5*(z/r)**2)) - 1/2*CD*S/m*rho*v_a*(xdot + thetadot*y);
    dX[4] = -mu*y*(1/r**3 + 3/2*J2*R**2/r**5*(1 - 5*(z/r)**2)) - 1/2*CD*S/m*rho*v_a*(ydot - thetadot*x);
    dX[5] = -mu*z*(1/r**3 + 3/2*J2*R**2/r**5*(3 - 5*(z/r)**2)) - 1/2*CD*S/m*rho*v_a*zdot
    dX[18:] = (A*phi).reshape((18*18,1))
    return dX

def Amatrix(X):
    R = 6378136.3
    S = 3
    m = 970
    rho0 = 3.614e-13
    r0 = 700000.0 + R
    H = 88667
    thetadot = 7.2921158543e-5

    x = X[0]
    y = X[1]
    z = X[2]
    xdot = X[3]
    ydot = X[4]
    zdot = X[5]
    xs = X[si[0]]
    ys = X[si[1]]
    zs = X[si[2]]

    r = sqrt(x**2 + y**2 + z**2)
    rho = rho0*np.exp(-(r - r0)/H)
    v_a = norm([xdot,ydot,zdot] - cross([0,0,thetadot],[x,y,z]))

    A = np.zeros((18,18))
    A[0][3] = 1
    A[1][4] = 1
    A[2][5] = 1
    A[3][0] = -mu/r**3*(1 - 3/2*J2*(R/r)**2*(5*(z/r)**2 - 1)) + 3*mu*x**2/r**5*(1 - 5/2*J2*(R/r)**2*(7*(z/r)**2 - 1)) + 1/2*CD*S/m*rho*(v_a*x*(xdot + thetadot*y)/(r*H) - (-thetadot*ydot + thetadot**2*x)*(xdot + thetadot*y)/v_a);
    A[3][1] = 3*mu*x*y/r**5*(1 - 5/2*J2*(R/r)**2*(7*(z/r)**2 - 1)) + 1/2*CD*S/m*rho*(v_a*y*(xdot + thetadot*y)/(r*H) - (thetadot*xdot + thetadot**2*y)*(xdot + thetadot*y)/v_a - v_a*thetadot);
    A[3][2] = 3*mu*x*z/r**5*(1 - 5/2*J2*(R/r)**2*(7*(z/r)**2 - 3)) + 1/2*CD*S/m*rho*v_a*z*(xdot + thetadot*y)/(r*H);
    A[3][3] = -1/2*CD*S/m*rho*((xdot + thetadot*y)**2/v_a + v_a);
    A[3][4] = -1/2*CD*S/m*rho*(ydot - thetadot*x)*(xdot + thetadot*y)/v_a;
    A[3][5] = -1/2*CD*S/m*rho*zdot*(xdot + thetadot*y)/v_a;
    A[3][6] = -x/r**3*(1 - 3/2*J2*(R/r)**2*(5*(z/r)**2 - 1));
    A[3][7] = 3/2*mu*x/r**3*(R/r)**2*(5*(z/r)**2 - 1);
    A[3][8] = -1/2*S/m*rho*v_a*(xdot + thetadot*y);
    A[4][0] = 3*mu*x*y/r**5*(1 - 5/2*J2*(R/r)**2*(7*(z/r)**2 - 1)) + 1/2*CD*S/m*rho*(v_a*(ydot - thetadot*x)/(r*H) - (thetadot**2*x - thetadot*ydot)*(ydot - thetadot*x)/v_a + v_a*thetadot);
    A[4][1] = -mu/r**3*(1 - 3/2*J2*(R/r)**2*(5*(z/r)**2 - 1)) + 3*mu*y**2/r**5*(1 - 5/2*J2*(R/r)**2*(7*(z/r)**2 - 1)) + 1/2*CD*S/m*rho*(v_a*y*(ydot - thetadot*x)/(r*H) - (thetadot*xdot + thetadot**2*y)*(ydot - thetadot*x)/v_a);
    A[4][2] = 3*mu*y*z/r**5*(1 - 5/2*J2*(R/r)**2*(7*(z/r)**2 - 3)) + 1/2*CD*S/m*rho*v_a*z*(ydot - thetadot*x)/(r*H);
    A[4][3] = -1/2*CD*S/m*rho*(ydot - thetadot*x)*(xdot + thetadot*y)/v_a;
    A[4][4] = -1/2*CD*S/m*rho*((ydot - thetadot*x)**2/v_a + v_a);
    A[4][5] = -1/2*CD*S/m*rho*zdot*(ydot - thetadot*x)/v_a;
    A[4][6] = -y/r**3*(1 - 3/2*J2*(R/r)**2*(5*(z/r)**2 - 1));
    A[4][7] = 3/2*mu*y/r**3*(R/r)**2*(5*(z/r)**2 - 1);
    A[4][8] = -1/2*S/m*rho*v_a*(ydot - thetadot*x);
    A[5][0] = 3*mu*x*z/r**5*(1 - 5/2*J2*(R/r)**2*(7*(z/r)**2 - 3)) + 1/2*CD*S/m*rho*(v_a*zdot*x/(r*H) - zdot*(thetadot**2*x - thetadot*ydot)/v_a);
    A[5][1] = 3*mu*y*z/r**5*(1 - 5/2*J2*(R/r)**2*(7*(z/r)**2 - 3)) + 1/2*CD*S/m*rho*(v_a*zdot*y/(r*H) - zdot*(thetadot*xdot + thetadot**2*y)/v_a);
    A[5][2] = -mu/r**3*(1 - 3/2*J2*(R/r)**2*(5*(z/r)**2 - 3)) + 3*mu*z**2/r**5*(1 - 5/2*J2*(R/r)**2*(7*(z/r)**2 - 5)) + 1/2*CD*S/m*rho*v_a*z*zdot/(r*H);
    A[5][3] = -1/2*CD*S/m*rho*zdot*(xdot + thetadot*y)/v_a;
    A[5][4] = -1/2*CD*S/m*rho*zdot*(ydot - thetadot*x)/v_a;
    A[5][5] = -1/2*CD*S/m*rho*(zdot**2/v_a + v_a);
    A[5][6] = -z/r**3*(1 - 3/2*J2*(R/r)**2*(5*(z/r)**2 - 3));
    A[5][7] = 3/2*mu*z/r**3*(R/r)**2*(5*(z/r)**2 - 3);
    A[5][8] = -1/2*S/m*rho*v_a*zdot;
    return A

def H_tilde(X,rho,rho_dot,theta,theta_dot,si):
    x = X[0]
    y = X[1]
    z = X[2]
    xdot = X[3]
    ydot = X[4]
    zdot = X[5]
    xs = X[si[0]]
    ys = X[si[1]]
    zs = X[si[2]]

    H_tilde = np.zeros(2,18)
    H_tilde[0][0] = (x - xs*np.cos(theta) + ys*np.sin(theta))/rho
    H_tilde[0][1] = (y - ys*np.cos(theta) - xs*np.sin(theta))/rho
    H_tilde[0][2] = (z - zs)/rho
    H_tilde[0][si[0]] = (xs - x*np.cos(theta) - y*np.sin(theta))/rho
    H_tilde[0][si[1]] = (ys - y*np.cos(theta) + x*np.sin(theta))/rho
    H_tilde[0][si[2]] = (zs - z)/rho
    H_tilde[1][0] = (xdot + thetadot*xs*np.sin(theta) + thetadot*ys*np.cos(theta))/rho - rhodot*(x - xs*np.cos(theta) + ys*np.sin(theta))/rho**2
    H_tilde[1][1] = (ydot + thetadot*ys*np.sin(theta) - thetadot*xs*np.cos(theta))/rho - rhodot*(y - ys*np.cos(theta) - xs*np.sin(theta))/rho**2
    H_tilde[1][2] = zdot/rho - rhodot*(z - zs)/rho**2
    H_tilde[1][3] = (x - xs*np.cos(theta) + ys*np.sin(theta))/rho
    H_tilde[1][4] = (y - ys*np.cos(theta) - xs*np.sin(theta))/rho
    H_tilde[1][5] = (z - zs)/rho
    H_tilde[1][si[0]] = (-xdot*np.cos(theta) + thetadot*x*np.sin(theta) - ydot*np.sin(theta) - thetadot*y*np.cos(theta))/rho - rhodot*(xs - x*np.cos(theta) - y*np.sin(theta))/rho**2
    H_tilde[1][si[1]] = (-ydot*np.cos(theta) + thetadot*y*np.sin(theta) + xdot*np.sin(theta) + thetadot*x*np.cos(theta))/rho - rhodot*(ys - y*np.cos(theta) + x*np.sin(theta))/rho**2
    H_tilde[1][si[2]] = -zdot/rho - rhodot*(zs - z)/rho**2
    return H_tilde

def kalman_filter(Xstar0, Y, x0bar, P0, R):
    thetadot = 7.2921158543e-5 #rad/s

    xhat = x0bar
    P = P0
    epsilon = np.zeros((len(Y),2))
    sol = propogator.propogate(ode,(0,24*3600),Xstar0)
    X = np.transpose(sol.y)
    
    for i in range(len(Y)):
        t = Y[i][0]
        station = Y[i][1]
        rho_obs = Y[i][2]
        rhodot_obs = Y[i][3]

        if station == 101:
            i_xs = 9
            i_ys = 10
            i_zs = 11
        elif station == 337:
            i_xs = 12
            i_ys = 13
            i_zs = 14
        else:
            i_xs = 15
            i_ys = 16
            i_zs = 17

        x = X[i][0]
        y = X[i][1]
        z = X[i][2]
        xdot = X[i][3]
        ydot = X[i][4]
        zdot = X[i][5]
        mu = X[i][6]
        J2 = X[i][7]
        CD = X[i][8]
        xs = X[i][i_xs]
        ys = X[i][i_ys]
        zs = X[i][i_zs]
        phi = reshape(X[i][18:],18,18)

        theta = thetadot*t
        rho = np.sqrt(x**2 + y**2 + z**2 + xz**2 + ys**2 + zs**2 - 2*(x*xs + y*ys)*np.cos(theta) + 2*(x*ys - y*xs)*np.sin(theta) - 2*z*zs)
        rhodot = np.sqrt(x*xdot + y*ydot + z*zdot - (xdot*xs + ydot*ys)*np.cos(theta) + thetadot*(x*xs + y*ys)*np.sin(theta) + (xdot*ys - ydot*xs)*np.sin(theta) + thetadot*(x*ys - y*xs)*np.cos(theta) - xdot*zs)/rho

        y = [rho_obs, rhodot_obs] - [rho, rhodot]
        H = H_tilde()
        K = P@np.transpose(H)@np.inverse(H@P@np.transpose(H) + R)

        xhat = xbar + K@(y - H@xbar)
        P = (np.eye(18) - K@H)@P@np.transpose(np.eye(18) - K@H) + K@R@np.transpose(K)
        Xhat = np.transpose(X[i][1:18}) + xhat
        epsilon[i] = np.transpose(y)

    sigma = std_dev(P)
    return sigma
