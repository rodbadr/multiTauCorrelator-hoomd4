"""
Test script for the multiTauCorrelator plugin in HOOMD 4.
This script sets up a simple simulation, applies the multiTauCorrelator action,
and verifies the output by comparing it with a C++ implementation.

It will compile the C++ codes in this directory and run them to compare the results
with the HOOMD implementation. 

The c++ codes are the same for both, this only checks if the hoomd implementation is working correctly
"""

import hoomd
import hoomd.md
import hoomd.logging
import unittest
import os
import logging
import errno
import subprocess
import numpy as np
from hoomd.multiTauCorrelator.analyze import autocorrelate
from hoomd.custom import Action
from hoomd.logging import log

# Compile Correlator_IO if not present or out of date
def compile_correlator_io():

    print("Compiling Correlator_IO...")
    ret = subprocess.call([
        "g++", "-O2", "-std=c++11", "-o", "Correlator_IO",
        "main_likh.cc", "test_correlator_likh.cc"
    ])
    if ret != 0:
        raise RuntimeError("Failed to compile Correlator_IO")

QUANTITIES = ["pressure_xy"]
FILENAME = "corr_otf_test.txt"

# method to separate the pressure tensor into its components and make them available as loggable quantities
class separatePressureTensor(hoomd.custom.Action):

    flags = [Action.Flags.PRESSURE_TENSOR]

    @log
    def pressure_xx(self):
        return self._state._simulation.operations.computes[0].pressure_tensor[0]

    @log
    def pressure_xy(self):
        return self._state._simulation.operations.computes[0].pressure_tensor[1]

    @log
    def pressure_xz(self):
        return self._state._simulation.operations.computes[0].pressure_tensor[2]

    @log
    def pressure_yy(self):
        return self._state._simulation.operations.computes[0].pressure_tensor[3]

    @log
    def pressure_yz(self):
        return self._state._simulation.operations.computes[0].pressure_tensor[4]

    @log
    def pressure_zz(self):
        return self._state._simulation.operations.computes[0].pressure_tensor[5]

    # this method is required for the Action class, but does nothing
    def act(self, timestep):
        pass


def silent_remove(filenames, disable=False):
    """
    Removes the target file name, catching and ignoring errors that indicate that the
    file does not exist.

    @param filename: The file to remove.
    @param disable: boolean to flag if want to disable removal
    """
    if not disable:
        for filename in filenames:
            try:
                os.remove(filename)
            except OSError as e:
                if e.errno != errno.ENOENT:
                    raise


# unit tests for autocorrelator

class test_correlate(unittest.TestCase):
    def setUp(self):
        self.device = hoomd.device.CPU()
        self.sim = hoomd.Simulation(device=self.device, seed=1)
        # Create a 2D square lattice manually using numpy
        n = [2, 2]  # 2x2 lattice
        a = 2.0
        N = n[0] * n[1]
        box_Lx = n[0] * a
        box_Ly = n[1] * a
        box_Lz = a
        positions = []
        for i in range(n[0]):
            for j in range(n[1]):
                positions.append([i * a + a/2 - box_Lx/2, j * a + a/2 - box_Ly/2, 0])
        positions = np.array(positions, dtype=np.float32)
        snap = hoomd.Snapshot()
        snap.particles.N = N
        snap.particles.position[:] = positions
        snap.particles.types = ['A']
        snap.particles.typeid[:] = 0
        snap.configuration.box = [box_Lx, box_Ly, box_Lz, 0, 0, 0]

        self.sim.create_state_from_snapshot(snap)

        self.thermo = hoomd.md.compute.ThermodynamicQuantities(hoomd.filter.All())
        self.sim.operations.computes.append(self.thermo)

        self.presTensor = separatePressureTensor()

        pres_writer = hoomd.write.CustomWriter(action=self.presTensor, trigger=hoomd.trigger.Periodic(1))
        self.sim.operations.writers.append(pres_writer)


    def test_create(self):
        """tests creation of the correlator """
        logger = hoomd.logging.Logger(categories=['scalar'])
        logger.add(self.presTensor, quantities=QUANTITIES[0])
        corr = autocorrelate(
            filename=FILENAME, quantities=QUANTITIES, logger=logger, eval_period=0
        )

    def test_create_file(self):
        """tests that the correct output file is created"""
        logger = hoomd.logging.Logger(categories=['scalar'])
        logger.add(self.presTensor, quantities=QUANTITIES[0])

        integrator = hoomd.md.Integrator(dt=0.01)
        langevin = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=1.0)
        integrator.methods.append(langevin)
        self.sim.operations.integrator = integrator
        corr = autocorrelate(
            filename=FILENAME, quantities=QUANTITIES, logger=logger, eval_period=0
        )
        correlator_writer = hoomd.write.CustomWriter(action=corr, trigger=hoomd.trigger.Periodic(1))
        self.sim.operations.writers.append(correlator_writer)
        self.sim.run(10)
        corr.write_to_file(self.sim.timestep)
        self.assertTrue(os.path.isfile(FILENAME))

    def test_values(self):
        """
        tests that the values in the output file are correct
        by comparing the hoomd implementation to the pure c++ implementation
        """
        logger = hoomd.logging.Logger(categories=['scalar'])
        logger.add(self.presTensor, quantities=QUANTITIES[0])
        
        integrator = hoomd.md.Integrator(dt=0.01)
        langevin = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=1.0)
        integrator.methods.append(langevin)
        self.sim.operations.integrator = integrator
        corr = autocorrelate(
            filename=FILENAME, quantities=QUANTITIES, logger=logger, eval_period=0
        )
        correlator_writer = hoomd.write.CustomWriter(action=corr, trigger=hoomd.trigger.Periodic(1))
        self.sim.operations.writers.append(correlator_writer)
        # Set up a logger to file for pressure
        file_logger = hoomd.logging.Logger(categories=['scalar', 'scalar'])
        file_logger.add(self.sim, quantities='timestep')
        file_logger.add(self.presTensor, quantities=QUANTITIES[0])
        
        with open("pressure_xy.log", "w", newline='\n') as f:
            table_file = hoomd.write.Table(output=f, trigger=hoomd.trigger.Periodic(1), logger=file_logger, delimiter=',', pretty=False, max_header_len=7)
            self.sim.operations.writers.append(table_file)
            self.sim.run(100)

        corr.write_to_file(self.sim.timestep)
        pressure_data = np.loadtxt("pressure_xy.log", skiprows=1, usecols=[1], delimiter=',')
        np.savetxt("pressure_data.txt", pressure_data)


        compile_correlator_io()

        subprocess.call(["./Correlator_IO", "pressure_data.txt", "corr_post_proc.txt"])

        corr_post_proc = np.loadtxt("corr_post_proc.txt")
        corr_otf = np.loadtxt("corr_otf_test.txt", skiprows=2, delimiter=",")
        np.testing.assert_almost_equal(corr_post_proc, corr_otf, decimal=2)
        silent_remove(["corr_post_proc.txt", "pressure_data.txt", "corr_otf_test.txt", "pressure_xy.log"])


if __name__ == "__main__":
    unittest.main(argv=["test_correlate.py", "-v"])

