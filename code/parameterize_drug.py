from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import ForceField
from pkg_resources import iter_entry_points

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, Draw, rdFMCS

import numpy as np

from simtk.openmm.app import PDBFile
from simtk.openmm import XmlSerializer

"""
Parameterizes the ES0186 using openforcefield for simulation with CB1

Also finds the maximum common substructure in the existing drug (MDMB-Fubinaca),
generates multiple conformations of ES0186 and aligns them all, choosing 
the best possible alignment as the starting coordinates.

"""

for entry_point in iter_entry_points(group='openforcefield.smirnoff_forcefield_directory'):
     print(entry_point.load()())
        
ff = ForceField('openff_unconstrained-1.0.0.offxml')

#New drug:
es0186 = Chem.MolFromSmiles('CCCCN1C=C(C(=O)N[C@H](C(N)=O)C(C)(C)C)C2=CC=CC=C12')
es0186=Chem.AddHs(es0186)

#Old drug:
mdmbfub_pdb = Chem.MolFromPDBFile('./processed_data/system_setup/MDMBFubinaca.pdb')
mdmbfub = Chem.MolFromSmiles('O=C(N[C@H](C(OC)=O)C(C)(C)C)C1=NN(CC2=CC=C(F)C=C2)C3=C1C=CC=C3')
mdmbfub_pdb = AllChem.AssignBondOrdersFromTemplate(mdmbfub, mdmbfub_pdb)

#Generate lots of conformers for new drug:
for _ in range(100):
    AllChem.EmbedMolecule(es0186, clearConfs=False)

##Find maximum common substructure so we can align based on this:
mcs =rdFMCS.FindMCS([es0186, mdmbfub_pdb])
smarts = mcs.smartsString
match = Chem.MolFromSmarts(smarts)
test_match_atoms = es0186.GetSubstructMatch(match)
ref_match_atoms = mdmbfub_pdb.GetSubstructMatch(match)

#Find alignments of all conformers of new drug to old drug:
alignments =[rdMolAlign.AlignMol(es0186, 
                    mdmbfub_pdb, 
                    prbCid=i,
                    atomMap=[[i,j] for i,j in zip(test_match_atoms, ref_match_atoms)]) for i in range(100)]

#Save best alignment:
Chem.MolToPDBFile(es0186, './processed_data/system_setup/drug.pdb', confId=int(np.argmin(alignments)))

##Now over to OpenForceField:
drug_pdbfile = PDBFile('./processed_data/system_setup/drug.pdb')
drug_mol = Molecule.from_smiles(Chem.MolToSmiles(es0186)) 
off_topology = Topology.from_openmm(openmm_topology=drug_pdbfile.topology,
                                   unique_molecules=[drug_mol])

#actual parameterizing step:
drug_system = ff.create_openmm_system(off_topology)

#Save as a serialized system:
with open('./processed_data/system_setup/drug_system.xml', 'w') as f:
    f.write(
            XmlSerializer.serialize(
                drug_system
            )
    )
