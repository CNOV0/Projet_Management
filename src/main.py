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
	CA = []
	with open(file_pdb, "r") as file_in:
		for ligne in file_in:
			if "ATOM" in ligne and "CA" in ligne:
				CA.append(ligne)
	data_ca = []
	for i in range(len(CA)):
		# coupe liste[i] selon les espaces et recup le n eme mot que l'on converti
		data = {"nb_ca" : int(CA[i].split()[5]), "x_ca": int(CA[i].split()[6]), \
				"y_ca": int(CA[i].split()[7]), "z_ca": int(CA[i].split()[8])}
		data_ca.append(data)
	return data_ca


def calculate_COM():
	# le vecteur position du COM d'un ensemble de solide vaut
	# la moyenne pondérée par les valeurs des masses des vect positions des solides de l'ensemble
	# Le COM : COM = (sum(miri) / M)
	# ou M est la somme des masses de particules; N est le nb de masses m;
	# ri est la distance du point de ref
	return 0


if __name__=="__main__":
	if len(sys.argv) != 2:
		sys.exit("Erreur : Il faut donner un seul argument qui est le nom du fichier pdb.")
	
	if not os.path.exists(sys.argv[1]):
		sys.exit("Erreur : Le fichier donné n'existe pas.")
	
	filename = sys.argv[1]
	data_dssp = dssp(filename)

