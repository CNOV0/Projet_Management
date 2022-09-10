import sys
import os
import math
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP


def dssp(file_pdb):
    """Run DSSP on a PDB file.

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
    """Get the value of accessibility for a given atom.

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


def find_ca_access(dssp_data, file_pdb):
    """Find the alpha carbon accessible to the solvant.

    Find the alpha carbon in a PDB file
    and select only those who are accessible to the solvant.

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
    carbon_alpha = []

    with open(file_pdb, "r") as file_in:
        for ligne in file_in:
            if ligne.startswith("ATOM") and "CA" in ligne:
                carbon_alpha.append(ligne)
            elif ligne.startswith("ENDMDL"):
                break

    ca_access = []

    for i in range(len(carbon_alpha)):
        access = access_solvant(dssp_data, i)
        if access > 0.5:
            data = [access, carbon_alpha[i].split()[3], \
                    float(carbon_alpha[i].split()[6]), \
                    float(carbon_alpha[i].split()[7]), \
                    float(carbon_alpha[i].split()[8])]
            ca_access.append(data)

    return pd.DataFrame(ca_access, columns=["accessibility", "resid_name",
                                            "x", "y", "z"])


def calculate_com(ca_data):
    """Calculate the center of mass of a protein.

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
    x_com = ca_data["x"].mean()
    y_com = ca_data["y"].mean()
    z_com = ca_data["z"].mean()

    return [x_com, y_com, z_com]


def center_protein(com, ca_data):
    """Center a protein at (0, 0, 0).

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
    ca_data["x"] -= com[0]
    ca_data["y"] -= com[1]
    ca_data["z"] -= com[2]


def fibonacci_sphere(samples=1000):
    """Generate uniforme distribution of points on a hemisphere.

    This function was taken from the internet.

    Parameters
    ----------
    samples : int
        The number of points wanted on the hemisphere

    Returns
    -------
    list of list
        The coordinates of each point on the hemisphere
    """
    points = []
    phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        coord_y = 1 - (i / float(samples - 1))  # y goes from 1 to 0
        radius = math.sqrt(1 - coord_y * coord_y)  # radius at y

        theta = phi * i  # golden angle increment

        coord_x = math.cos(theta) * radius
        coord_z = math.sin(theta) * radius

        points.append([coord_x, coord_y, coord_z])

    return points


def calculate_hydrobicity(sphere_points, ca_data):
    """Calculate the hydrophobicity within a slice.

    Parameters
    ----------
    sphere_points : list of list
        The coordinates of points on a hemisphere
    ca_data : Pandas DataFrame
        The coordinates of alpha carbon of a protein and their residu's name

    The DataFrame also contains the accessibility of each atom.

    Returns
    -------
    int
        The best score of hydrophobity
    """
    score = 0
    mini = pd.Series([ca_data["x"].min(), ca_data["y"].min(), ca_data["z"].min()])
    maxi = pd.Series([ca_data["x"].max(), ca_data["y"].max(), ca_data["z"].max()])
    # Loop on each point of the sphere
    for point in sphere_points:
        # should be parallelised
        while mini["x"] < maxi["x"] or mini["y"] < maxi["y"] or mini["z"] < maxi["z"]:
            for i in range(ca_data.shape[0]):
                nb_ca, hydro = 0, 0
                if carbon_is_in_membrane(point, ca_data.iloc[i], mini, 22):
                    nb_ca += 1
                    if is_hydrophobe(ca_data.iloc[i, "resid_name"]):
                        hydro += 1
                score_slice = hydro / nb_ca
                if score < score_slice:
                    score = score_slice  # add info to locate it
            mini += 1
    return score


def is_hydrophobe(resid_name):
    """Test if an atom is in a hydrophobe residu.

    Parameters
    ----------
    resid_name : str
        The name of the atom's residu

    Returns
    -------
    Boolean
        If the atom is or not in a hydrophobe residu.
    """
    hydrophobe = ["TRP", "ISO", "LEU", "PHE", "ALA", "MET", "VAL"]

    if resid_name in hydrophobe:
        return True
    return False


def carbon_is_in_membrane(sphere_pt, ca_coords, point_min, memb_width):
    """Test if an atom is between two planes.

    Parameters
    ----------
    sphere_pt : list
        The coordinate of a point on the hemisphere
    ca_coords : Pandas Series
        The coordinate of an atom
    point_min : Pandas Series
        The coordinate of the lowest point
    memb_width : int
        The width of the membrane

    Returns
    -------
    Boolean
        If the atom is or not between the planes.
    """
    plane1_d = calculate_plane(sphere_pt, point_min, 0)
    plane2_d = calculate_plane(sphere_pt, point_min, memb_width)
    ca_d = calculate_plane(sphere_pt, ca_coords, 0)

    if plane1_d < ca_d < plane2_d:
        return True
    return False


def calculate_plane(sphere_pt, atom, memb_width):
    """Calculate the plane perpendicular to the sphere point's axis.

    Parameters
    ----------
    sphere_point : list
        The coordinate of a point on the hemisphere
    atom : Pandas Series
        The coordinate of an atom
    memb_width : int
        The width of the membrane

    Returns
    -------
    list of tuple
        The coordinates of each point on the hemisphere
    """
    plane_d = sphere_pt[0] * (atom["x"] + sphere_pt[0] * memb_width) + \
            sphere_pt[1] * (atom["y"] + sphere_pt[1] * memb_width) + \
            sphere_pt[2] * (atom["z"] + sphere_pt[2] * memb_width)
    return plane_d


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Erreur : Il faut donner un seul argument qui est"
                " le nom du fichier pdb.")

    if not os.path.exists(sys.argv[1]):
        sys.exit("Erreur : Le fichier donné n'existe pas.")

    filename = sys.argv[1]

    dssp_info = dssp(filename)
    ca_info = find_ca_access(dssp_info, filename)
    center_of_mass = calculate_com(ca_info)
    center_protein(center_of_mass, ca_info)
    uniform_points = fibonacci_sphere(10)
