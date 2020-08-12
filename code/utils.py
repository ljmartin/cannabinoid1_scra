from rdkit import Chem
import pandas as pd
from rdkit.Chem import PandasTools

def get_molecules():
   df = pd.read_csv('./data/smiles.csv')
   PandasTools.AddMoleculeColumnToFrame(df,'smiles','molecule',includeFingerprints=False)
   return df
