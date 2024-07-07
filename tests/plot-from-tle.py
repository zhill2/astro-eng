import context
import numpy as np
from tools import *
from planets import earth

tle = [
    "WORLDVIEW-2 (WV-2)",
    "1 35946U 09055A   24064.55783087  .00000225  00000+0  89652-4 0  9997",
    "2 35946  98.4212 139.4686 0000796 200.0130 160.1026 14.37503040755752"
]
state = tle_importer.tle2state(tle, earth.mu)
    
lambda_model = lambda t, X : models.base_model(earth.mu, t, X)
    
start = 0
end = 5*24*3600

ode_sol = propogator.propogate(lambda_model, (start, end), state)
states = np.array(ode_sol.y)
plotter.plot_orbit(states, surf=plotter.create_surf(earth.R))
