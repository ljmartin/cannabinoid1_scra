import parmed
import numpy as np
import itertools

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

from sys import stdout
import sys
"""
Makes a complex with the CHARMM-GUI system and our newly parameterized drug. 

Using parmed because then it's easy to just 'add' structures together in python.
i.e. structure3 = structure1 + structure2
"""

ID = str(sys.argv[1])

#this is required to recognise DPPC molecules:
topology.Topology.loadBondDefinitions('./data/dppc.xml')

#Load up CHARMM-GUI pdb, and create a PARMED structure from it.  
prot_pdb = PDBFile('./data/step5_assembly.pdb')
omm_forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml', 'amber14/lipid17.xml')
prot_system = omm_forcefield.createSystem(prot_pdb.topology, rigidWater=False)
prot_structure = parmed.openmm.load_topology(prot_pdb.topology,
                                           prot_system,
                                           xyz=prot_pdb.positions)


#Load up the parameterized drug system, and again make it into a parmed structure:
drug_system = XmlSerializer.deserialize(open('./processed_data/off_systems/drug_system_'+ID+'.xml').read())
drug_pdbfile = PDBFile('./processed_data/aligned_drugs/drug_'+ID+'.pdb')

drug_structure = parmed.openmm.load_topology(drug_pdbfile.topology,
                                                drug_system,
                                                xyz=drug_pdbfile.positions)

#This is the biggest step but it takes 1 second:
complex_structure = prot_structure + drug_structure

#Now give this new structure some periodic box conditions:
#first, shift so that minimum coordinate is at (0,0,0). This allows to set PBC around it.
complex_structure.coordinates = complex_structure.coordinates-np.min(complex_structure.coordinates, axis=0)
#set PBC with a short buffer:
newbox = (np.max(complex_structure.coordinates,axis=0)+np.array(3))*angstroms
complex_structure.box = (newbox[0], newbox[0], newbox[2], 90, 90, 90)
#now return coordinates so that protein COM is at (0,0,0). This is (unfortunately)
#required to make the barostat play nicely with the restraints during equilibration.
protein_indices = np.array([count for count, i in enumerate(prot_pdb.topology.atoms()) if i.residue.name not in ["DPPC", "CL", "NA", "HOH"]])
complex_structure.coordinates = complex_structure.coordinates-np.mean(complex_structure.coordinates[protein_indices], axis=0)

#Turn into an OpenMM System object for simulations:
#These settings will be stuck unless you re-run this script! Luckily they're pretty standard settings.
complex_system = complex_structure.createSystem(nonbondedMethod=PME,
                                                nonbondedCutoff=0.9*nanometer,
                                               constraints=HBonds,
                                                rigidWater=True)

#Save output:
complex_structure.save('./processed_data/complex_systems/complex_coords_'+ID+'.pdb', overwrite=True)

#PSF files don't like having numbered atom types, because at 10,000 VMD fails
#So just set them all to zero.
for a in complex_structure.atoms:
    a.type = '0'
complex_structure.save('./processed_data/complex_systems/complex_struct_'+ID+'.psf', overwrite=True)

with open('./processed_data/complex_systems/complex_system_'+ID+'.xml', 'w') as f:
    f.write(
            XmlSerializer.serialize(
                complex_system
            )
    )
