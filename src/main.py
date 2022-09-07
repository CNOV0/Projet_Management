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


def access_solvant():
	# recup tous les acces Puis selectionner que ceux qui sont sup a 1 seuil AVEC INDEX
	a_key = list(data_dssp.keys())[1]
	data_dssp[a_key]


def find_CA(file_pdb):
	# trouve CA Puis selection de ceux qui sont en accord avec acces select
	with open(file_pdb, "r") as file_in:
		for ligne in file_in:
			piou = 0


if __name__=="__main__":
	if len(sys.argv) != 2:
		sys.exit("Erreur : Il faut donner un seul argument qui est le nom du fichier pdb.")
	
	if not os.path.exists(sys.argv[1]):
		sys.exit("Erreur : Le fichier donné n'existe pas.")
	
	filename = sys.argv[1]
	data_dssp = dssp(filename)

