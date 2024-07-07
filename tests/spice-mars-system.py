import matplotlib.pyplot as plt
import numpy as np
import spiceypy as spice
import urllib.request

urllib.request.urlretrieve('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls', './spice/latest_leapseconds.tls')
urllib.request.urlretrieve('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp', './spice/de440.bsp')
urllib.request.urlretrieve('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc','./spice/gm_de440.tpc')
urllib.request.urlretrieve('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc','./spice/pck00011.tpc')

spice.furnsh('./spice/latest_leapseconds.tls')
spice.furnsh('./spice/de440.bsp')
spice.furnsh('./spice/gm_de440.tpc')
spice.furnsh('./spice/pck00011.tpc')

class Body:
    def __init__(self,ID,currpos=[0, 0, 0]):
        self.ID = ID
        self.name = spice.bodc2n(int(ID))
        self.radius = spice.gdpool("BODY"+ID+"_RADII",0,1)[0]
        self.min_altitude = 0
        self.GM = spice.gdpool("BODY"+ID+"_GM",0,1)[0]
        self.J2 = None
        self.currpos = currpos
        self.create_surf()

        print("Creating Object", ID)
        print("Name:",self.name)
        print("Radius:",self.radius)
        print("GM:",self.GM)
        print("Current Position:",self.currpos)
        print()

    def create_surf(self):
        pi = np.pi
        cos = np.cos
        sin = np.sin
        phi, theta = np.mgrid[0.0:pi:20j, 0.0:2.0*pi:20j]
        x = self.radius*sin(phi)*cos(theta)
        y = self.radius*sin(phi)*sin(theta)
        z = self.radius*cos(phi)

        x += self.currpos[0]
        y += self.currpos[1]
        z += self.currpos[2]

        self.surf = [x, y, z]

    def get_position(self, et, cb, ref='J2000', abcorr='LT+S'):
        pos, t = spice.spkpos(self.ID, et, ref, abcorr, cb.ID)
        return pos.reshape((3,1))
    
    def orbit(self, start, dur, obs, ref='J2000', abcorr='LT+S'):
        print("Running", self.name, "orbit for", dur, "seconds centered at",obs.name)
        if type(start) == str:
            start = spice.str2et(start)
        tspan = np.linspace(start,start+dur,num=1000)
        self.positions, self.t = spice.spkpos(self.ID, tspan, ref, abcorr, obs.ID)
        self.currpos = self.positions[-1]
        self.positions = self.positions.T
        print("Updated position to", self.currpos)
        self.create_surf()

if __name__=='__main__':
    urllib.request.urlretrieve('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/mar097.bsp','./spice/mar097.bsp')
    spice.furnsh('./spice/mar097.bsp')

    mars = Body("499")
    phobos = Body("401")
    deimos = Body("402")

    start = '2024 02 29'
    dur = 3600*24*10

    phobos.orbit(start, dur, mars)
    deimos.orbit(start, dur, mars)

    fig = plt.figure(figsize=(9, 9))
    ax  = fig.add_subplot(111, projection='3d')
    ax.plot(phobos.positions[0], phobos.positions[1], phobos.positions[2])
    ax.plot(deimos.positions[0], deimos.positions[1], deimos.positions[2])
    ax.plot_surface(mars.surf[0], mars.surf[1], mars.surf[2])
    ax.plot_surface(phobos.surf[0], phobos.surf[1], phobos.surf[2])
    ax.plot_surface(deimos.surf[0], deimos.surf[1], deimos.surf[2])
    plt.title('Mars System')
    plt.show()
