import os
import sys
import math
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

def dssp(file_pdb):
	"""Run DSSP on a PDB file

	Used to get the solvant accessibility of each atom.

	Parameters
	----------
	file_pdb : str
		The name of a PDB file

		The filename may contain de path to fetch the file.

	Returns
	-------
	list ?
		Different informations on the atom
		- especially the accessibility.
	"""
    parsed = PDBParser()
    structure = parsed.get_structure("2K1A", file_pdb)
    model = structure[0]
    dssp_data = DSSP(model, file_pdb)
    
    return dssp_data


def access_solvant(dssp_data, index):
	"""Get the value of accessibility for a given atom

	Parameters
	----------
	dssp_data : list ?
		A list containing various informations and the accessibility
	index : int
		The atom for which we want the accessibility

	Returns
	-------
	int
		The solvant accessibility of the atom
	"""
	a_key = list(dssp_data.keys())[index]
	return dssp_data[a_key][3]


def find_CA_access(dssp_data, file_pdb):
	"""Find the alpha carbon accessible to the solvant

	Find the alpha carbon in a PDB file
	and select only those who are accessible to the solvant

	Parameters
	----------
	dssp_data : list ?
		A list containing various informations and the accessibility
	file_pdb : str
		The name of a PDB file

		The filename may contain de path to fetch the file.

	Returns
	-------
	Pandas DataFrame
		DataFrame containing the solvant accessibility, the residu's name,
		the coordinates of the alpha carbon.
	"""
	CA = []
	
	with open(file_pdb, "r") as file_in:
		for ligne in file_in:
			if ligne.startswith("ATOM") and "CA" in ligne:
				CA.append(ligne)
			elif ligne.startswith("ENDMDL"):
				break
			
	ca_access = []
	
	for i in range(len(CA)):
		access = access_solvant(dssp_data, i)
		if access > 0.5:
			data = [access, CA[i].split()[3], \
					float(CA[i].split()[6]), \
					float(CA[i].split()[7]), \
					float(CA[i].split()[8])]
			ca_access.append(data)
	
	return pd.DataFrame(ca_access, columns=["accessibility", "resid_name", "x", "y", "z"])


def calculate_COM(ca_data):
	"""Calculate the center of mass of a protein

	Parameters
	----------
	ca_data : Pandas DataFrame
		DataFrame containing the atoming coordinates of an alpha carbon
        
		The DataFrame contains also other informations such as
		the solvant accessibility and the residu's name.

	Returns
	-------
	list
        List containing the center of mass of a
        protein for each atomic coordinate
	"""
	x_COM = ca_data["x"].mean()
	y_COM = ca_data["y"].mean()
	z_COM = ca_data["z"].mean()
	
	return [x_COM, y_COM, z_COM]


def center_protein(COM, ca_data):
	"""Center a protein in (0, 0, 0)

	Parameters
	----------
	COM: list
		List containing the atomic coordinates of
		the center of mass of a protein
	ca_data : Pandas DataFrame
		DataFrame containing the coordinates of the alpha carbon of a protein

		The DataFrame contains also other information such as
		the solvant accessibility of the atom and the residu's name.

	Returns
	-------
	int
		Le produit des deux nombres.
	"""
	ca_data["x"] -= COM[0]
	ca_data["y"] -= COM[1]
	ca_data["z"] -= COM[2]


def fibonacci_sphere(samples=1000):
	"""Generate uniforme distribution of points on a hemisphere

	This function was taken from the internet.

	Parameters
	----------
	samples : int
		The number of points wanted on the hemisphere

	Returns
	-------
	list of tuple
		The coordinates of each point on the hemisphere
	"""
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


def calculate_hydrobicity(sphere_points, ca_coords):
	membrane_width = 22  # For protein with just CA - 15 A otherwise
	point_min = [ca_coords["x"].min(), ca_coords["y"].min(), ca_coords["z"].min()]
	point_max = [ca_coords["x"].max(), ca_coords["y"].max(), ca_coords["z"].max()]
	# Loop on each point of the sphere
	for point in sphere_points:
		# should be parallelised
		while point_min < point_max:
			for ligne in ca_coords:  # pas sure
				nb_ca, hydro = 0, 0
				if carbon_is_in_membrane(point, ARG):
					nb_ca += 1
					if ca_coords.loc[INDICE, "resid_name"] is in [HYDROPHOBES]:
						hydro += 1
			hydro /= nb_ca
			#  keep best one w/ info to locate it
			point_min = 0  # new mean point (add 1A)
	return hydrophobicity


def carbon_is_in_membrane():
	plane1_d = calcule_plane()
	plane2_d = calculate_plane()
	ca_d = calculate_plane()  # sphere_points[0] * atom.x + sphere_points[1] * atom.y + sphere_points[2] * atom.z
	if (ca_d > plane1_d) and (ca_d < plane2_d):
		return True
	else:
		return False


def calculate_plane():
	x_0 = ca_data.loc[0, "x"]
	# calcule l'équation de plan
	plane_d = sphere_points[0] * (min_point[0] + sphere_points[0] * membraine_width) + \
			sphere_points[1] * (min_point[1] + sphere_point[1] * membraine_width) + \
			sphere_points[2] * (min_point[2] + sphere_points[2] * membraine_width)
	return plane_d


if __name__=="__main__":
	if len(sys.argv) != 2 :
		sys.exit("Erreur : Il faut donner un seul argument qui est le nom du fichier pdb.")
	
	if not os.path.exists(sys.argv[1]):
		sys.exit("Erreur : Le fichier donné n'existe pas.")
	
	filename = sys.argv[1]
	
	dssp_data = dssp(filename)
	ca_access = find_CA_access(data_dssp, filename)
	COM = calculate_COM(ca_access)
	center_protein(COM, ca_access)
	sphere_points = fibonacci_sphere(10)

