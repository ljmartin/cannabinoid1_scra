from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import ForceField
from pkg_resources import iter_entry_points

from rdkit import Chem

import utils

from simtk.openmm.app import PDBFile
from simtk.openmm import XmlSerializer

for entry_point in iter_entry_points(group='openforcefield.smirnoff_forcefield_directory'):
     print(entry_point.load()())

ff = ForceField('openff_unconstrained-1.0.0.offxml')


def make_off_system(mol, ID):##Now over to OpenForceField:
    drug_pdbfile = PDBFile('./processed_data/aligned_drugs/drug_'+str(ID)+'.pdb')
    drug_mol = Molecule.from_smiles(Chem.MolToSmiles(mol))
    off_topology = Topology.from_openmm(openmm_topology=drug_pdbfile.topology,
                                   unique_molecules=[drug_mol])

    #actual parameterizing step:
    drug_system = ff.create_openmm_system(off_topology)
    
    return drug_system

def serialize(ds, ID):
    #Save as a serialized system:
    with open('./processed_data/off_systems/drug_system_'+str(ID)+'.xml', 'w') as f:
        f.write(
            XmlSerializer.serialize(
                ds
            )
        )


if __name__=="__main__":
    df = utils.get_molecules()
    for count, row in df.iterrows():
        mol = row['molecule']
        ID = row['id']

        drug_system = make_off_system(mol, ID)
        serialize(drug_system, ID)
