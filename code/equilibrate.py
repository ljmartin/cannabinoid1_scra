from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

import numpy as np

from sys import stdout
import sys
from tqdm import tqdm



"""
Equilibrates the drug-protein complex. Starts with strong restraints on both protein 
and drug and gradually releases them. The drug restraint is released
a little bit faster so it can find it's niche before the protein moves to 
accommodate it's starting structure. 
"""

ID = sys.argv[1]

#load
complex_pdb = PDBFile('./processed_data/complex_systems/complex_coords_'+ID+'.pdb')
complex_system = XmlSerializer.deserialize(open('./processed_data/complex_systems/complex_system_'+ID+'.xml').read())

#get prot and drug indices for doing restraints later. 
protein_indices = np.array([count for count, i in enumerate(complex_pdb.topology.atoms()) if i.residue.name not in ["DPPC", "DPP", "CLA", "NA", "HOH"]])
drug_indices = np.array([count for count, i in enumerate(complex_pdb.topology.atoms()) if i.residue.name=='UNL'])

##add a barostat. Not only will this initially fix any vacuum between periodic cells,
##it will be maintained throughout to keep pressure reasonable. 
barostat = MonteCarloMembraneBarostat(1*bar, 200*bar*nanometer, 310*kelvin, MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree)
complex_system.addForce(barostat)


#A restraining force on the drug atoms for equilibration:
proteinForce = CustomExternalForce("proteinF*periodicdistance(x, y, z, x0, y0, z0)^2")
proteinForce.addGlobalParameter("proteinF", 5*kilocalories_per_mole/angstroms**2)
proteinForce.addPerParticleParameter("x0")
proteinForce.addPerParticleParameter("y0")
proteinForce.addPerParticleParameter("z0")
for count, coord in enumerate(complex_pdb.positions):
    if count in protein_indices:
        proteinForce.addParticle(int(count), [coord[0], coord[1], coord[2]])

print("Does it use periodic conditions?", proteinForce.usesPeriodicBoundaryConditions())
complex_system.addForce(proteinForce)


##A restraining force on the drug atoms for equilibration:
#drugForce = CustomExternalForce("drugF*periodicdistance(x, y, z, x0, y0, z0)^2")
#drugForce.addGlobalParameter("drugF", 5*kilocalories_per_mole/angstrom**2)
#drugForce.addPerParticleParameter("x0")
#drugForce.addPerParticleParameter("y0")
#drugForce.addPerParticleParameter("z0")
#for count, coord in enumerate(complex_pdb.positions):
#    if count in drug_indices:
#        drugForce.addParticle(int(count), [coord[0], coord[1], coord[2]])
#
#print("Does it use periodic conditions?", drugForce.usesPeriodicBoundaryConditions())
#complex_system.addForce(drugForce)


#OpenMM simulation machinery:
integrator = LangevinMiddleIntegrator(310*kelvin, 10/picosecond, 0.003*picoseconds)
platform = Platform.getPlatformByName('OpenCL')
prop = {'OpenCLPrecision':'single', 'OpenCLDeviceIndex':'0'}
simulation = Simulation(complex_pdb.topology, complex_system, integrator, platform, prop)

#minimize starting coords:
simulation.context.setPositions(complex_pdb.positions)
simulation.minimizeEnergy(maxIterations=150)
#for temp in [150, 175, 200, 225, 250, 275, 300, 310]:
#    integrator.setTemperature(150*kelvin)
#    simulation.step(1000)
    
#add a DCD writer:
simulation.reporters.append(DCDReporter('./processed_data/equilibration/traj_equilibration_'+ID+'.dcd', 10000))
simulation.reporters.append(StateDataReporter('./processed_data/equilibration/log_equilibration_'+ID+'.log', 1000, 
                                              step=True, 
                                              potentialEnergy=True, 
                                              temperature=True, speed=True))


#reduce restraints:
#we are looking to avoid clashes more than remove any structure errors
#from crystallization, since receptor activity changes take order 10^3-4 nanoseconds
for i in tqdm(range(300)):
    simulation.step(4000) #12ps per it, = 12*300, 3ns.
    #simulation.context.setParameter('drugF', simulation.context.getParameter('drugF')*0.9)
    simulation.context.setParameter('proteinF', simulation.context.getParameter('proteinF')*0.9)

##The restraint forces are basically zero now, so remove both:
forces = np.array([isinstance(i, CustomExternalForce) for i in complex_system.getForces()])
while np.any(forces):
    complex_system.removeForce(int(forces.nonzero()[0][0]))
    forces = np.array([isinstance(i, CustomExternalForce) for i in complex_system.getForces()])

##Extra 1 ns for good measure.
simulation.step(400000)



##Save it all.
#System:
with open('./processed_data/equilibration/equilibration_system_'+ID+'.xml', 'w') as f:
    f.write(
            XmlSerializer.serialize(
                complex_system
            )
    )

#Current state:
with open('./processed_data/equilibration/equilibration_state_'+ID+'.xml', 'w') as f:
    f.write(
            XmlSerializer.serialize(
                simulation.context.getState(getPositions=True,
                                     getForces=True, getEnergy=True,
                                     enforcePeriodicBox=True)
            )
    )

#PDB file (we already have one but this just keeps the three in the equilibration dir for neatness):
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('./processed_data/equilibration/equilibration_coords_'+ID+'.pdb', 'w'))
