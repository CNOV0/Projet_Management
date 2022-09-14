# Transmembrane Protein Project

Author : Maya Zygadlo

Here is the code to find the most probable localisation of the membrane from a Protein Data Bank file.

## Setup your environment

Clone the repository:

```bash
git clone https://github.com/CNOV0/Projet_Management.git
```

Move to the new directory:

```bash
cd Projet_Management
```

Install [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Install [mamba](https://github.com/mamba-org/mamba):

```bash
conda install mamba -n base -c conda-forge
```

Create the "Projet_TMP" conda environment:
```
mamba env create -f src/Projet_TMP.yaml
```
or if you don't want to install mamba:
```
conda env create -f src/Projet_TMP.yaml
```

Load the "Projet_TMP" conda environment:
```
conda activate Projet_TMP
```

Note: you can also update the conda environment with:

```bash
mamba env update -f Projet_TMP.yaml
```

To deactivate an active environment, use

```
conda deactivate
```

## Run the algorithm

Create the results directory:
```bash
mkdir results
```

Run the algorithm:
```bash
python3 src/main.py [file]
```
