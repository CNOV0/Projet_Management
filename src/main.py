from Bio.PDB.PDBParser import PDBParser
par = PDBParser()
structure = par.get_structure('2K1A','Projet_Membrane_Maya/data/2k1a.pdb')
import Bio.PDB.NACCESS as nac  # module utiliser pour trouver les molecules accessibles au solvant

tmp_dir = "Projet_Membrane_Maya/results/tmp/"
Bio.PDB.NACCESS.run_naccess(model=structure[0], \
                            pdb_file="Projet_Membrane_Maya/data/2k1a.pdb", \
                            temp_path=tmp_dir)
