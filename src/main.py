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


def access_solvant(data_dssp, index):
	a_key = list(data_dssp.keys())[index]
	return data_dssp[a_key][3]


def find_CA_access(data_dssp, file_pdb):
	# trouve CA Puis selection de ceux qui sont en accord avec acces select
	CA = []
	with open(file_pdb, "r") as file_in:
		for ligne in file_in:
			if ligne.startswith("ATOM") and "CA" in ligne:
				CA.append(ligne)
			elif ligne.startswith("ENDMDL"):
				break
			
	ca_access = []
	print(len(CA))
	for i in range(len(CA)):
		access = access_solvant(data_dssp, i)
		if access > 0.5:
			data = {"access" : access, "nom_resid" : CA[i].split()[3], "x_ca": float(CA[i].split()[6]), \
					"y_ca": float(CA[i].split()[7]), "z_ca": float(CA[i].split()[8])}
			ca_access.append(data)
	return ca_access


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
	ca_access = find_CA_access(data_dssp, filename)
	print(ca_access)

