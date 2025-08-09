# On-the-fly autocorrelation for HOOMD-blue version 4

## About

The multiTauCorrelator Custom Action for HOOMD-blue v4 provides on-the-fly time autocorrelation of any HOOMD quantity using the [Likhtman multiple tau autocorrelation algorithm](https://aip.scitation.org/doi/10.1063/1.3491098).

## Installation

To use this Custom Action, clone the repository into the directory of choice then build and install it as follows:

    git clone https://github.com/atravitz/correlator-hoomd-plugin
    cd correlator-hoomd-plugin
    cmake -B build/correlator_plugin/
    cmake --build build/correlator_plugin/
    cmake --install build/correlator_plugin/

Afterwards, it is available as a HOOMD4 module `hoomd.multiTauCorrelator`.

The Custom Action can be tested by running the included test script:

    python3 test/test_correlator.py

## Documentation


multiTauCorrelator is a HOOMD4 Custom Action, as such, users should have a basic understanding of creating and running HOOMD scripts.

    class hoomd.multiTauCorrelator.analyze.autocorrelate(
        quantities : list[str],
        logger : hoomd.logging.Logger,
        filename : str="autocorrelate.log",
        eval_period : int=0,
        numcorrin : int=32,
        p_in : int=16,
        m_in : int=2,
        normalize : bool=False)

Parameters:
* `quantities : list[str]` - List of quantities to autocorrelate (must be present in the logger)
* `logger : hoomd.logging.Logger` - Logger object containing the quantities to be autocorrelated
* `filename : str` - File to write autocorrelation values to
* `eval_period : int` - Autocorrelation data is written every eval_period time steps
* `numcorrin : int` - Number of correlators (levels) in the multi-tau correlator hierarchy
* `p_in : int` - Number of time bins for each correlator
* `m_in : int` - Decimation factor (averaging factor) between correlator levels
* `normalize : bool` - If True, normalizes the correlation function by subtracting the average value at 'tau=0' from all values

Additional autocorrelate functions are:

    autocorrelate.write_to_file(timestep)

Writes out autocorrelation values to file at the given timestep

Quantities that can be autocorrelated are any built-in or user-defined quantities that can be logged in HOOMD-blue.

A full example is provided in the `examples/example_correlator.py` file. In the example we calculate custom quantities and use them with the autocorrelate Custom Action to compute their time correlation functions.

A quick example:

```python

import hoomd
import hoomd.md
import hoomd.logging
import numpy as np
from hoomd.multiTauCorrelator.analyze import autocorrelate

# Set up the simulation device and create a Simulation object
device = hoomd.device.CPU()
sim = hoomd.Simulation(device=device, seed=1)

# Create a 2D square lattice of particles using numpy
n = [2, 2]  # 2x2 lattice
a = 2.0  # lattice spacing
N = n[0] * n[1]
box_Lx = n[0] * a
box_Ly = n[1] * a
box_Lz = a
positions = []

# Generate lattice positions centered in the box
for i in range(n[0]):
    for j in range(n[1]):
        positions.append([i * a + a/2 - box_Lx/2, j * a + a/2 - box_Ly/2, 0])

positions = np.array(positions, dtype=np.float32)

# Create a HOOMD snapshot and initialize particle positions and box
snap = hoomd.Snapshot()
snap.particles.N = N
snap.particles.position[:] = positions
snap.particles.types = ['A']
snap.particles.typeid[:] = 0
snap.configuration.box = [box_Lx, box_Ly, box_Lz, 0, 0, 0]

sim.create_state_from_snapshot(snap)

# Add a compute for thermodynamic quantities (e.g., pressure)
thermo = hoomd.md.compute.ThermodynamicQuantities(hoomd.filter.All())
sim.operations.computes.append(thermo)

# Set up a logger to record the pressure
logger = hoomd.logging.Logger(categories=['scalar'])
logger.add(thermo, quantities=['pressure'])

# Set up the integrator and Langevin thermostat
integrator = hoomd.md.Integrator(dt=0.01)
langevin = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=1.0)
integrator.methods.append(langevin)
sim.operations.integrator = integrator

# Create the autocorrelate Custom Action object for pressure
corr = autocorrelate(
    filename='correlate.log', quantities=['pressure'], logger=logger
)

# Add the correlator to the simulation writers so it runs every step
correlator_writer = hoomd.write.CustomWriter(action=corr, trigger=hoomd.trigger.Periodic(1))
sim.operations.writers.append(correlator_writer)

# Run the simulation for 100 timesteps
sim.run(100)

# Write the correlation data to file
corr.write_to_file(sim.timestep)

```

## Acknowledgement

- The original C++ code was obtained from:

    https://blogs.upm.es/compsoftmatter/software/multiple-tau-correlator/

    based on the algorithm published in:

    Ramírez, J., Sukumaran, S.K., Vorselaars, B. and Likhtman, A.E., 2010. Efficient on the fly calculation of time correlation functions in computer simulations. J. Chem. Phys. 133, 154103 (2010)

    DOI: 10.1063/1.3491098



- Inspired by the plugin for HOOMD-blue v2 by Alyssa Travitz and Alexander Adams:

    https://github.com/atravitz/correlator-hoomd-plugin

- HOOMD4 Custom Action by Rodrique Badr