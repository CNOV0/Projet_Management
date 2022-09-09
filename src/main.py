import os
import sys
import math
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
	CA = []
	
	with open(file_pdb, "r") as file_in:
		for ligne in file_in:
			if ligne.startswith("ATOM") and "CA" in ligne:
				CA.append(ligne)
			elif ligne.startswith("ENDMDL"):
				break
			
	ca_access = []
	
	for i in range(len(CA)):
		access = access_solvant(data_dssp, i)
		if access > 0.5:
			data = {"access" : access, "nom_resid" : CA[i].split()[3], \
					"x_ca": float(CA[i].split()[6]), \
					"y_ca": float(CA[i].split()[7]), \
					"z_ca": float(CA[i].split()[8])}
			ca_access.append(data)
	
	return ca_access


def calculate_COM(ca_data):
	x_COM, y_COM, z_COM = 0, 0, 0
	
	for ca in ca_data:
		x_COM += ca["x_ca"]
		y_COM += ca["y_ca"]
		z_COM += ca["z_ca"]
		
	x_COM /= len(ca_data)
	y_COM /= len(ca_data)
	z_COM /= len(ca_data)
	
	return [x_COM, y_COM, z_COM]


def center_protein(COM, ca_data):
	for ca in ca_data:
		ca["x_ca"] -= COM[0]
		ca["y_ca"] -= COM[1]
		ca["z_ca"] -= COM[2]


def fibonacci_sphere(samples=2000):
	points = []
	phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians

	for i in range(samples):
		y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
		radius = math.sqrt(1 - y * y)  # radius at y

		theta = phi * i  # golden angle increment

		x = math.cos(theta) * radius
		z = math.sin(theta) * radius

		if y >= 0:
			points.append([x, y, z])

	return points


if __name__=="__main__":
	if len(sys.argv) != 2 :
		sys.exit("Erreur : Il faut donner un seul argument qui est le nom du fichier pdb.")
	
	if not os.path.exists(sys.argv[1]):
		sys.exit("Erreur : Le fichier donné n'existe pas.")
	
	filename = sys.argv[1]
	
	data_dssp = dssp(filename)
	ca_access = find_CA_access(data_dssp, filename)
	COM = calculate_COM(ca_access)
	center_protein(COM, ca_access)
	uniform_points = fibonacci_sphere(20)

