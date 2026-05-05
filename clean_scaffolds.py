#read data

import sqlite3
import pandas as pd

# Connect to the database
conn = sqlite3.connect('/DISCO1/Aylin/Data_Bases/ChEMBL/chembl_36/chembl_36.db')

# Read the SQL query results into a Pandas DataFrame
df = pd.read_sql_query('SELECT * FROM compound_structures', conn)

# Close the connection
conn.close()

from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem

def extract_scaffold_and_gf(s):
    mol = Chem.MolFromSmiles(s)
    if mol is None:
        return (None, None)
    
    # 1. Strip Salts (keeps the largest fragment)
    remover = SaltRemover()
    mol = remover.StripMol(mol, dontRemoveEverything=True)
    
    # 2. Neutralize Charges
    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)
    
    # 3. Canonicalize Tautomers
 #   te = rdMolStandardize.TautomerEnumerator()
 #   mol = te.Canonicalize(mol)
    
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    scaffold_smiles = Chem.MolToSmiles(scaffold) if scaffold is not None else None
    
    try:
        fw = MurckoScaffold.MakeScaffoldGeneric(scaffold) if scaffold is not None else None
        gf_smiles = Chem.MolToSmiles(fw) if fw is not None else None
    except Exception:
        gf_smiles = None
    return (scaffold_smiles, gf_smiles)


# Apply the function once and expand the result into two columns
# Note: this uses pandas.Series to turn the tuple into separate columns
df[['scaffold', 'generic_framework']] = df['canonical_smiles'].apply(lambda s: pd.Series(extract_scaffold_and_gf(s)))

print(df.head())

df[['molregno','canonical_smiles', 'scaffold', 'generic_framework']].to_csv('~/Chemical_Space/chembl_scaffolds.csv', index=False)
