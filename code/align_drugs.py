import utils
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, Draw, rdFMCS
import numpy as np



def generate_conformers(mol, n=100):
    """Generates 3d structures for a molecule using Riniker's ETKDG method. N=number of conformers to make"""
    try:
        for _ in range(n):
            AllChem.EmbedMolecule(mol, clearConfs=False)
    except:
        print('Failed to embed')



def align_mcs(new_mol, ref_mol):
    """Ref mol is the crystallized ligand. New mol is the drug you want to add to the structure.
    It returns an array of scores (one for each conformer), where the lowest score is best"""
    ##Find maximum common substructure so we can align based on this:
    mcs =rdFMCS.FindMCS([new_mol, ref_mol])
    smarts = mcs.smartsString
    match = Chem.MolFromSmarts(smarts)
    test_match_atoms = new_mol.GetSubstructMatch(match)
    ref_match_atoms = ref_mol.GetSubstructMatch(match)
    
    #Find alignments of all conformers of new drug to old drug:
    alignments_scores =[rdMolAlign.AlignMol(new_mol,
                    ref_mol,
                    prbCid=i,
                    atomMap=[[i,j] for i,j in zip(test_match_atoms, ref_match_atoms)]) for i in range(100)]
    
    return alignments_scores



if __name__=="__main__":
    #old drug in crystal structure is MDMB Fubinaca:
    mdmbfub_pdb = Chem.MolFromPDBFile('./data/MDMBFubinaca.pdb')
    mdmbfub = Chem.MolFromSmiles('O=C(N[C@H](C(OC)=O)C(C)(C)C)C1=NN(CC2=CC=C(F)C=C2)C3=C1C=CC=C3')
    mdmbfub_pdb = AllChem.AssignBondOrdersFromTemplate(mdmbfub, mdmbfub_pdb)

    #get new drugs:
    df = utils.get_molecules()
    #for each drug, align to the co-crystallized ligand, and save PBD:
    for count, row in df.iterrows():
        idx = row['id']
        mol = row['molecule']
        mol = Chem.AddHs(mol)
        generate_conformers(mol, n=100)
    
        alignments = align_mcs(mol, mdmbfub_pdb)
    
        Chem.MolToPDBFile(mol, './processed_data/aligned_drugs/drug_'+str(idx)+'.pdb', confId=int(np.argmin(alignments)))
