"""
Example script to run a simulation with DPD forces and calculate the
autocorrelations required for stress relaxation function calculations.

The output file will be named `correlation.csv` and will contain the
autocorrelation data for the specified quantities.

The output file will have the following format:
```

correlator evaluated at timestep {timestep}
timestep, corr_pxy, corr_pxz, corr_pyz, corr_Nxy, corr_Nxz, corr_Nyz
t0,corr_pxy1,corr_pxz1,...
t1,corr_pxy2,corr_pxz2,...
...

```
"""

import hoomd
import hoomd.md
import numpy as np
from hoomd.multiTauCorrelator.analyze import autocorrelate

from hoomd.custom import Action
from hoomd.logging import log

pi = np.pi


class separatePressureTensor(hoomd.custom.Action):

    """
    Custom Action to separate the pressure tensor into its components.
    The difference between the diagonal components Nab is also calculated.
    This is useful for stress relaxation function calculations.
    """

    flags = [Action.Flags.PRESSURE_TENSOR]

    @log
    def pxx(self):
        return self._state._simulation.operations.computes[0].pressure_tensor[0]

    @log
    def pxy(self):
        return self._state._simulation.operations.computes[0].pressure_tensor[1]

    @log
    def pxz(self):
        return self._state._simulation.operations.computes[0].pressure_tensor[2]

    @log
    def pyy(self):
        return self._state._simulation.operations.computes[0].pressure_tensor[3]

    @log
    def pyz(self):
        return self._state._simulation.operations.computes[0].pressure_tensor[4]

    @log
    def pzz(self):
        return self._state._simulation.operations.computes[0].pressure_tensor[5]
    
    @log
    def Nxy(self):
        pxx = self._state._simulation.operations.computes[0].pressure_tensor[0]
        pyy = self._state._simulation.operations.computes[0].pressure_tensor[3]
        return pxx - pyy
    
    @log
    def Nxz(self):
        pxx = self._state._simulation.operations.computes[0].pressure_tensor[0]
        pzz = self._state._simulation.operations.computes[0].pressure_tensor[5]
        return pxx - pzz

    @log
    def Nyz(self):
        pyy = self._state._simulation.operations.computes[0].pressure_tensor[3]
        pzz = self._state._simulation.operations.computes[0].pressure_tensor[5]
        return pyy - pzz


    def act(self, timestep):
        # This method is called at each timestep, but we don't need to do anything here
        # so we just pass
        pass





##################################################
### Set parameters and read initial conditions ###
##################################################

buf = 0.075 # buffer for neighbour list

T = 1.0; # temperature for Langevin Thermostat
rc = 1.0; # cutoff radius for dpd attraction
gamma = 4.5 # drag coefficient
dt = 1e-2 # time step for the integrator

A = 40 # DPD force parameter


nsteps = 1000 # number of steps to run the simulation

device = hoomd.device.CPU() # Create a CPU device for the simulation
sim = hoomd.Simulation(device=device, seed=1) # Create a simulation object with a random seed
# Create a 2D square lattice manually using numpy

# 3D cubic lattice parameters
n = [2, 2, 2]  # 2x2x2 cubic lattice
a = 2.0 # lattice spacing
N = n[0] * n[1] * n[2] # total number of particles
box_Lx = n[0] * a # box length in x direction
box_Ly = n[1] * a # box length in y direction
box_Lz = n[2] * a # box length in z direction
positions = []
# Generate positions for the cubic lattice
# not very efficient, but works for small lattices
for i in range(n[0]):
    for j in range(n[1]):
        for k in range(n[2]):
            positions.append([
                i * a + a/2 - box_Lx/2,
                j * a + a/2 - box_Ly/2,
                k * a + a/2 - box_Lz/2
            ])

# Convert positions to numpy array
positions = np.array(positions, dtype=np.float32)
snap = hoomd.Snapshot() # Create a snapshot to initialize the simulation state
snap.particles.N = N # Set the number of particles in the snapshot
snap.particles.position[:] = positions # Set the particle positions in the snapshot
snap.particles.types = ['A'] # Set the availabe particle types
snap.particles.typeid[:] = 0 # Set the particle type IDs (all particles are of type 'A')
snap.configuration.box = [box_Lx, box_Ly, box_Lz, 0, 0, 0] # Set the simulation box dimensions in the snapshot

sim.create_state_from_snapshot(snap) # Initialize the simulation state from the snapshot

#neighbour list via cell
nl = hoomd.md.nlist.Cell(buffer = buf, exclusions = ())

#disipative particle dynamics force
dpd = hoomd.md.pair.DPD(nlist=nl, kT=T, default_r_cut=rc)
dpd.params[('A', 'A')] = dict(A=A, gamma=gamma) # Set the DPD force parameters for the particle type 'A'

all = hoomd.filter.All() # Create a filter to select all particles

integrator = hoomd.md.Integrator(dt=dt) # Create an integrator with a time step of dt

integrator.forces.append(dpd) # Add the DPD force to the integrator

nve = hoomd.md.methods.ConstantVolume(filter=all) # Create a constant volume method to integrate the equations of motion
integrator.methods.append(nve) # Append the constant volume method to the integrator

sim.operations.integrator = integrator # Set the integrator for the simulation

sim.state.thermalize_particle_momenta(filter=all, kT=T)

# Compute thermodynamic properties such as pressure and temperature
# create the compute object for thermodynamic quantities and append it to the simulation compute operations
thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=all)
sim.operations.computes.append(thermodynamic_properties)

# Create a custom action object to separate the pressure tensor into its components
# This will make all the loggers and methods available as attributes of the object
presTensor = separatePressureTensor()

# Create a custom writer for the pressure tensor components at each timestep
# although our logger is not writing anything, this will make sure that the 
# separated components are available at every timestep and can be used for
# further analysis or logging
pres_writer = hoomd.write.CustomWriter(action=presTensor, trigger=hoomd.trigger.Periodic(1))
sim.operations.writers.append(pres_writer) # Append the pressure writer to the simulation writers

# List of quantities to correlate
quantities = ['pxy', 'pxz', 'pyz', 'Nxy', 'Nxz', 'Nyz']

### Create a logger to log the quantities of interest
# first create a logger object with the desired categories
stressRel_logger = hoomd.logging.Logger(categories=["scalar"])
# add the pressure tensor components to the logger. This will make sure that the quantities are logged
stressRel_logger.add(presTensor, quantities=quantities)

# Create the autocorrelator action object with the specified quantities and logger
# this is the class we defined in the custom action file `analyze.py`
correlator_action = autocorrelate(quantities=quantities, logger=stressRel_logger, filename=f"correlation.csv")

# Create a custom writer for the correlator action.
# Although the action itself does not write anything,
# this will make sure that the act() method is called at each trigger timestep
correlator_writer = hoomd.write.CustomWriter(action=correlator_action, trigger=hoomd.trigger.Periodic(1))
sim.operations.writers.append(correlator_writer)

print('Simulation Starting')

sim.run(nsteps)

# call the write_to_file method from the correlator_action object
# this will write the correlation data to the specified file
correlator_action.write_to_file(sim.timestep)

print('Simulation ended')
