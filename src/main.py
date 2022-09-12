import sys
import os
import math
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import calculate_membrane
import extract_ca_accessible
import center


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Erreur : Il faut donner un seul argument qui est"
                " le nom du fichier pdb.")

    if not os.path.exists(sys.argv[1]):
        sys.exit("Erreur : Le fichier donné n'existe pas.")

    filename = sys.argv[1]

    dssp_info = extract_ca_accessible.dssp(filename)
    ca_info = extract_ca_accessible.find_ca_access(dssp_info, filename)
    center_of_mass = center.calculate_com(ca_info)
    center.center_protein(center_of_mass, ca_info)
    uniform_points = calculate_membrane.fibonacci_sphere()
    memb_pos = calculate_membrane.position(uniform_points, ca_info)
    calculate_membrane.plot(ca_info, memb_pos, center_of_mass)
    print("The membrane plane has been calculated.", 
    "\nYou can view the result in the results directory.")
