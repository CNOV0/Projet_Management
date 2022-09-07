from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
p = PDBParser()
structure = p.get_structure("2K1A", "Projet_Membrane_Maya/data/2k1a.pdb")
model = structure[0]
dssp = DSSP(model, "Projet_Membrane_Maya/data/2k1a.pdb")

a_key = list(dssp.keys())[1]
dssp[a_key]
