import os
import sys
import math
import pandas as pd
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
			data = [access, CA[i].split()[3], \
					float(CA[i].split()[6]), \
					float(CA[i].split()[7]), \
					float(CA[i].split()[8])]
			ca_access.append(data)
	
	return pd.DataFrame(ca_access, columns=["accessibility", "nom_resid", "x", "y", "z"])


def calculate_COM(ca_data):
	x_COM = ca_data["x"].mean()
	y_COM = ca_data["y"].mean()
	z_COM = ca_data["z"].mean()
	
	return [x_COM, y_COM, z_COM]


def center_protein(COM, ca_data):
	ca_data["x"] -= COM[0]
	ca_data["y"] -= COM[1]
	ca_data["z"] -= COM[2]


def fibonacci_sphere(samples=1000):
	points = []
	phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians

	for i in range(samples):
		y = 1 - (i / float(samples - 1))  # y goes from 1 to 0
		radius = math.sqrt(1 - y * y)  # radius at y

		theta = phi * i  # golden angle increment

		x = math.cos(theta) * radius
		z = math.sin(theta) * radius
		
		points.append((x, y, z))

	return points


def atom_is_transmembrane(sphere_points, ca_coords):
	membrane_width = 22  # Dans le cas d'une protein ayant just eun squelette CA
	point_min = [0, 0, 0]  # valeur a mettre (xmin,ymin,zmin)
	point_max = [0, 0, 0]
	# boucle sur chaque point sphere et check si est dans plane
	for point in sphere_points:
		# ici paraleliser - itere sur min -> max par saut de 1A
		while point_min < point_max:
			if carbon_is_in_membrane():
				hydro = 0
				# check si hydrophobe
				# calcul hydrophobicité
				# garder best one w/ info pour la retrouvé
			point_min = 0  # new mean point (ajout d'1A)
	return hydrophobicity


def carbon_is_in_membrane():
	d_plane1 = calcule_plane()
	d_plane2 = calculate_plane()
	d_ca = calculate_plane()  # sphere_points[0] * atom.x + sphere_points[1] * atom.y + sphere_points[2] * atom.z
	if (d_ca > d_plane1) and (d_ca < d_plane2):
		return True
	else:
		return False


def calculate_plane():
	x_0 = ca_data.loc[0, "x"]
	# calcule l'équation de plan
	d_plan = sphere_points[0] * (min_point[0] + sphere_points[0] * membraine_width) + \
			sphere_points[1] * (min_point[1] + sphere_point[1] * membraine_width) + \
			sphere_points[2] * (min_point[2] + sphere_points[2] * membraine_width)
	return d_plan


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
	sphere_points = fibonacci_sphere(10)

