import os
from netCDF4 import Dataset
import numpy as np
import opendrift.models.openoil.tamoc_chemical_properties as chem
# from tamoc import seawater


class Plume(object):
    """
    Loads the netCDF file with the results from the TAMOC bent plume model

    Translates the state space variables for a Lagrangian plume element into
    its individual parts and derived quantitites.

    Stores the results as plume object attributes and tamoc_dbm.FluidPArticle
    objects.

    Parameters
    ----------
    plume_file : string
        The name of the netCDF file holding the Tamoc bent plume simulation
        results.

    Attributes
    ----------
    ntracers : int
        Number of passive chemical tracers (--)
    nchems : int
        Number of chemicals tracked for dissolution of the dispersed phase
        particles (--)
    nparticles : int
        Number of dispersed phase particles (--)
    particles : dict of tamoc_dispersed_phases.TrackedFluidParticle objects
        For integer key: the Fluid Particle object holding the derived
        properties of the plume particles at the end of the near field
        simulation.
    t : float
        Independent variable for the current time (s)
    q : ndarray
        Dependent variable for the current state space
    M : float
        Mass of the Lagrangian element (kg)
    Se : float
        Salt in the Lagrangian element (psu kg)
    He : float
        Heat of the Lagrangian element (J)
    Jx : float
        Dynamic momentum of the Lagrangian element in the x-direction
        (kg m/s)
    Jy : float
        Dynamic momentum of the Lagrangian element in the y-direction
        (kg m/s)
    Jz : float
        Dynamic momentum of the Lagrangian element in the z-direction
        (kg m/s)
    H : float
        Relative thickness of the Lagrangian element h/V (s)
    x : float
        Current x-position of the Lagrangian element (m)
    y : float
        Current y-position of the Lagrangian element (m)
    z : float
        Current z-position of the Lagrangian element (m)
    s : float
        Current s-position along the centerline of the plume for the
        Lagrangian element (m)
    M_p : dict of ndarrays
        For integer key: the total mass fluxes (kg/s) of each component in a
        particle.
    H_p : ndarray
        Total heat flux for each particle (J/s)
    t_p : ndarray
        Time since release for each particle (s)
    X_p : ndarray
        Position of each particle in local plume coordinates (l,n,m) (m).
    cpe : ndarray
        Masses of the chemical components involved in dissolution (kg)
    cps : ndarray
        Masses of dissolved chemicals per second released at the end of the plume simulation (kg/s)
    cte : ndarray
        Masses of the passive tracers in the plume (concentration kg)
    Pa : float
        Ambient pressure at the current element location (Pa)
    Ta : float
        Ambient temperature at the current element location (K)
    Sa : float
        Ambient salinity at the current element location (psu)
    ua : float
        Crossflow velocity in the x-direction at the current element location
        (m/s)
    ca_chems : ndarray
        Ambient concentration of the chemical components involved in
        dissolution at the current element location (kg/m^3)
    ca_tracers :
        Ambient concentration of the passive tracers in the plume at the
        current element location (concentration)
    rho_a : float
        Ambient density at the current element location (kg/m^3)
    S : float
        Salinity of the Lagrangian element (psu)
    T : float
        Temperature of the Lagrangian element (T)
    c_chems :
        Concentration of the chemical components involved in dissolution for
        the Lagrangian element (kg/m^3)
    c_tracers :
        Concentration of the passive tracers in the Lagrangian element
        (concentration)
    u : float
        Velocity in the x-direction of the Lagrangian element (m/s)
    v : float
        Velocity in the y-direction of the Lagrangian element (m/s)
    w : float
        Velocity in the z-direction of the Lagrangian element (m/s)
    hvel : float
        Velocity in the horizontal plane for the Lagrangian element (m/s)
    V : float
        Velocity in the s-direction of the Lagrangian element (m/s)
    h : float
        Current thickness of the Lagrangian element (m)
    rho : float
        Density of the entrained seawater in the Lagrangian element (kg/m^3)
    b : float
        Half-width of the Lagrangian element (m)
    sin_p : float
        The sine of the angle phi (--)
    cos_p : float
        The cosine of the angle phi (--)
    sin_t : float
        The sine of the angle theta (--)
    cos_t : float
        The cosine of the angle theta (--)
    phi : float
        The vertical angle from horizontal of the current plume trajectory
        (rad in range +/- pi/2).  Since z is positive down (depth), phi =
        pi/2 point down and -pi/2 points up.
    theta : float
        The lateral angle in the horizontal plane from the x-axis to the
        current plume trajectory (rad in range 0 to 2 pi)
    mp : ndarray
        Masses of each component in a particle (Kg).

    """
    def __init__(self, chem_data_file=None, plume_file=None):
        super(Plume, self).__init__()
        # Create a Plume object from a saved file
        if chem_data_file:
            self.load_chem_data(chem_data_file)
        else:
            __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__), 'data'))
            chem_data = os.path.join(__location__, 'pseudo_chemdata.csv')
            self.load_chem_data(chem_data)
        if plume_file:
            self.plume_file = plume_file
            self.load_plume(plume_file)

    def load_chem_data(self, fname=None):
        self.chem_data, self.chem_units = chem.load_data(fname)
        self.chem_names = self.chem_data.keys()

    def load_plume(self, fname=None):
        # Open the netCDF dataset object containing the simulation results
        nc = Dataset(fname)

        required_dims = ('ns', 'nchems', 'nparticles')
        missing_dims = [d for d in required_dims if d not in nc.dimensions]
        if missing_dims:
            nc.close()
            raise ValueError(
                f'TAMOC plume file {fname} missing required dimensions: '
                f'{missing_dims}')

        required_vars = ('t', 'nb0', 'nbe', 'fp_type', 'q')
        missing_vars = [v for v in required_vars if v not in nc.variables]
        if missing_vars:
            nc.close()
            raise ValueError(
                f'TAMOC plume file {fname} missing required variables: '
                f'{missing_vars}')

        # Compute the dimensions of the arrayed data
        ns = len(nc.dimensions['ns'])
        self.nchems = len(nc.dimensions['nchems'])

        self.nparticles = len(nc.dimensions['nparticles'])
        nt = nc.variables['t'].n_times

        # Extract the arrayed data
        self.t = np.zeros(nt)
        self.t[:] = nc.variables['t'][0:nt,0]
        self.nb0 = np.zeros(self.nparticles)
        self.nb0[:] = nc.variables['nb0'][0:self.nparticles]
        self.nbe = np.zeros(self.nparticles)
        self.nbe[:] = nc.variables['nbe'][0:self.nparticles]
        self.fp_type = np.zeros(self.nparticles)
        self.fp_type[:] = nc.variables['fp_type'][0:self.nparticles]
        dt = self.nbe[0] / self.nb0[0]
        q = np.zeros((nt, ns))
        for i in range(ns):
            q[:,i] = nc.variables['q'][0:nt,i]

        # Extract the state-space variables from q
        self.q = q
        self.M = q[:,0]
        self.Se = q[:,1]
        self.He = q[:,2]
        self.Jx = q[:,3]
        self.Jy = q[:,4]
        self.Jz = q[:,5]
        self.H = q[:,6]
        self.x = q[:,7]
        self.y = q[:,8]
        self.z = q[:,9]
        idx = 11
        M_p = {}     #  np.zeros([self.nparticles, self.nchems], dtype=np.float32)
        M_p0 = {}     # np.zeros([self.nparticles, self.nchems], dtype=np.float32)
        self.mp = np.zeros([self.nparticles, self.nchems], dtype=np.float32)
        self.mp0 = np.zeros([self.nparticles, self.nchems], dtype=np.float32)
        for i in range(self.nparticles):
            m_c = np.zeros([self.nchems], dtype=np.float32)

            m_c0 = np.zeros([self.nchems], dtype=np.float32)
            for j in range(self.nchems):
                if q[-1, idx + j] > 0.:
                    m_c[j] = q[-1, idx + j]
                if q[0, idx + j] > 0.: m_c0[j] = q[0, idx + j]
            M_p[i] = m_c / dt
            M_p0[i] = m_c0 / dt
            idx += self.nchems
            self.mp[i, :] = M_p[i] / self.nb0[i]
            self.mp0[i, :] = M_p0[i] / self.nb0[i]
            idx += 5
        self.M_p = M_p
        self.cpe = q[-1, idx:idx + self.nchems]
        self.cps = self.cpe / dt
        idx += self.nchems
        if len(q[-1,idx:]) >= 1:
            self.ntracers = len(q[-1,idx:])
            self.cte = q[-1, idx:]
        else:
            self.ntracers = 0
            self.cte = np.array([])

        # Compute the derived quantities
        self.S = self.Se / self.M


    def get_particle_properties(self, particles_file, tracked=False):

        self.kg_s_ini = np.sum(self.mp0, axis=1) * self.nb0
        self.kg_s_diss_plume = np.sum(self.mp0 - self.mp, axis=1) * self.nb0
        self.kg_s_end = np.sum(self.mp, axis=1) * self.nb0
        self.tracked = tracked

        if tracked:
            t_max = np.zeros(self.nparticles)
            z_end_track = np.zeros(self.nparticles)
            mass_exit_plume = np.zeros((self.nparticles, self.nchems))
            mass_end_track = np.zeros((self.nparticles, self.nchems))
            mass_surface = np.zeros((self.nparticles, self.nchems))
            for i in range(self.nparticles):
                idx = '000'[:3-len(str(i))] + str(i)
                particle_nc = Dataset(self.plume_file[:-3] + idx + '.nc')
                t_max[i] = particle_nc.variables['t'][-1]
                z_end_track[i] = particle_nc.variables['y'][-1, 2]
                mass_exit_plume[i,:] = particle_nc.variables['y'][0, 3:-1]
                mass_end_track[i,:] = particle_nc.variables['y'][-1, 3:-1]
                if z_end_track[i] < 10.:
                    mass_surface[i,:] = mass_end_track[i,:]
            self.kg_s_end = np.sum(mass_end_track, axis=1) * self.nb0
            self.kg_s_diss_farfield = np.sum(mass_exit_plume - mass_end_track, axis=1) * self.nb0
            self.kg_s_surface = np.sum(mass_surface, axis=1) * self.nb0

        # Open the netCDF dataset object containing the particles' positions
        nc = Dataset(particles_file)
        # Extract the particle positions
        self.xp = np.zeros(self.nparticles)
        self.xp[:] = nc.variables['xp'][0:self.nparticles]
        self.yp = np.zeros(self.nparticles)
        self.yp[:] = nc.variables['yp'][0:self.nparticles]
        self.zp = np.zeros(self.nparticles)
        self.zp[:] = nc.variables['zp'][0:self.nparticles]
        self.tp = np.zeros(self.nparticles)
        self.tp[:] = nc.variables['tp'][0:self.nparticles]
        self.t_hyd = np.zeros(self.nparticles)
        self.t_hyd[:] = nc.variables['t_hyd'][0:self.nparticles]
        self.T0 = np.zeros(self.nparticles)
        self.T0[:] = nc.variables['T0'][0:self.nparticles]

        # Compute the final values of the derived quantities at the end of the simulation
        self.rho = 1035          # seawater.density(self.T[-1], self.S[-1], Pa)
        self.c_chems = self.cpe / (self.M[-1] / self.rho)
        self.c_tracers = self.cte / (self.M[-1] / self.rho)
        self.u = self.Jx[-1] / self.M[-1]
        self.v = self.Jy[-1] / self.M[-1]
        self.w = self.Jz[-1] / self.M[-1]
        self.hvel = np.sqrt(self.u ** 2 + self.v ** 2)
        self.V = np.sqrt(self.hvel ** 2 + self.w ** 2)
        self.h = self.H[-1] * self.V
        self.b = np.sqrt(self.M[-1] / (self.rho * np.pi * self.h))
        self.sin_p = self.w / self.V
        self.cos_p = self.hvel / self.V
        if self.hvel == 0.:
            # if hvel = 0, flow is purely along z; let theta = 0
            self.sin_t = 0.
            self.cos_t = 1.
        else:
            self.sin_t = self.v / self.hvel
            self.cos_t = self.u / self.hvel
        self.phi = np.arctan2(self.w, self.hvel)
        self.theta = np.arctan2(self.v, self.u)