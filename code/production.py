from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, exit, stderr
import numpy as np
from tqdm import tqdm



ID = sys.argv[1]

complex_pdb = PDBFile('./processed_data/equilibration/equilibration_coords_'+ID+'.pdb')

complex_system = XmlSerializer.deserialize(open('./processed_data/equilibration/equilibration_system_'+ID+'.xml').read())

for count, force in enumerate(complex_system.getForces()):
    print(force)


statefile = './processed_data/equilibration/equilibration_state_'+ID+'.xml'
state = XmlSerializer.deserialize(open(statefile).read())

integrator = LangevinMiddleIntegrator(323*kelvin, 1/picosecond, 0.003*picoseconds)
platform = Platform.getPlatformByName('OpenCL')
prop = {'OpenCLPrecision':'single', 'OpenCLDeviceIndex':'0'}
simulation = Simulation(complex_pdb.topology, complex_system, integrator, platform, prop)

simulation.reporters.append(DCDReporter('./processed_data/production/traj_'+ID+'.dcd', 50000))
simulation.reporters.append(StateDataReporter('./processed_data/production/log_'+ID+'.log', 10000,
                                              step=True,
                                              time=True,
                                              potentialEnergy=True,
                                              temperature=True,
                                              totalSteps=5000000,
                                              speed=True))


#simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
#simulation.context.setPositions(state.getPositions(asNumpy=True)+1*nanometers)
#simulation.context.setVelocitiesToTemperature(310*kelvin)
print('loadstate')
simulation.loadState(statefile)

simulation.step(30000000)
#for nanosecond in tqdm(range(30)):
#    simulation.step(100000)


#Current state:  
with open('./processed_data/production/state_'+ID+'.xml', 'w') as f:
    f.write(
            XmlSerializer.serialize(
                simulation.context.getState(getPositions=True,
                                     getForces=True, getEnergy=True,
                                     enforcePeriodicBox=True)
            )
    )
