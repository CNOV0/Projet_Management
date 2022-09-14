#!/usr/bin/python3

"Module to calculate the centre of mass and center a protein."

# import modules
import pandas as pd


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
    ca_data["x"] -= com[0]  # on every x of the dtf
    ca_data["y"] -= com[1]
    ca_data["z"] -= com[2]
