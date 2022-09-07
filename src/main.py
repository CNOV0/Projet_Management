import os
import sys
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

def dssp(file_pdb):
    parsed = PDBParser()
    structure = parsed.get_structure("2K1A", file_pdb)
    model = structure[0]
    data_dssp = DSSP(model, file_pdb)
    return data_dssp


if __name__=="__main__":
	if len(sys.argv) != 2:
		sys.exit("Erreur : Il faut donner un seul argument qui est le nom du fichier pdb.")

	if not os.path.exists(sys.argv[1]):
		sys.exit("Erreur : Le fichier donné n'existe pas.")
	
	filename = sys.argv[1]
	data_dssp = dssp(filename)
	a_key = list(data_dssp.keys())[1]
	print(data_dssp[a_key])

