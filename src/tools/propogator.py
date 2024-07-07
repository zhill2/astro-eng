from scipy.integrate import solve_ivp

def propogate(model, tspan, x0, method='RK45'):
    """
    Propogates the orbit for the given timespan. Wraps the scipy ivp integrator.

    Parameters
    ----------
    model  : The function for the ode. It must only have two inputs the current
             time step and the current state.
    tspan  : A tuple containing the start and end time
    x0     : The initial state of the system
    method : The method to use to solve the IVP. (Default is EK45)

    Returns
    -------
    Bunch object with the following fields defined:
    t : ndarray, shape (n_points,)
        Time points.
    y : ndarray, shape (n, n_points)
        Values of the solution at `t`.
    sol : `OdeSolution` or None
        Found solution as `OdeSolution` instance; None if `dense_output` was
        set to False.
    t_events : list of ndarray or None
        Contains for each event type a list of arrays at which an event of
        that type event was detected. None if `events` was None.
    y_events : list of ndarray or None
        For each value of `t_events`, the corresponding value of the solution.
        None if `events` was None.
    nfev : int
        Number of evaluations of the right-hand side.
    njev : int
        Number of evaluations of the Jacobian.
    nlu : int
        Number of LU decompositions.
    status : int
        Reason for algorithm termination:

            * -1: Integration step failed.
            *  0: The solver successfully reached the end of `tspan`.
            *  1: A termination event occurred.

    message : string
        Human-readable description of the termination reason.
    success : bool
        True if the solver reached the interval end or a termination event
        occurred (``status >= 0``).
    """
    ode_sol = solve_ivp(
        fun = model,
        t_span = tspan,
        y0 = x0,
        method = method,
        rtol = 1e-12,
        atol = 1e-12,
        dense_output = False
    )
    return ode_sol
