# Project 2 - Clustering and Drug Discovery
## Ryder Easterlin
## Due 02/19/2021

![BuildStatus](https://github.com/rle1323/Project2/workflows/HW2/badge.svg?event=push)

In this assignment, you will evaluate results from a high-throughput virtual screen against the SARS-CoV2 Spike protein / Human ACE2 interface.  There are two parts to this assignment and Part 2 requires completion of Part 1. We recommend reading through both Part 1 and Part 2 before beginning this assignment. 

* Part 1 - API and implementation
* Part 2 - Evaluating clustering

The data we are considering comes from [Smith and Smith, 2020](https://chemrxiv.org/articles/preprint/Repurposing_Therapeutics_for_the_Wuhan_Coronavirus_nCov-2019_Supercomputer-Based_Docking_to_the_Viral_S_Protein_and_Human_ACE2_Interface/11871402). In this study, they generated 6 Spike-Ace2 interface poses using MD simulations. They then docked ~10k small molecules against each protein conformation. Provided for you is the top (#1) pose for each ligand docked against one Spike-ACE2 interface conformation, as well as the corresponding SMILES string, AutoDock Vina score, and the “On” bits in the Extended Connectivity Fingerprint for that compound. These can all be found in ligand\_information.csv.


### main
Runs all code in align/\_\_main\_\_.py, useful for part 2
```
python -m algs
```

### testing
Testing is as simple as running
```
python -m pytest test/*
```
from the root directory of this project.

## Documentation

### Functions

read_ligand_data(filename)
```
Reads in the provided ligand data from csv format, and stores the ligands as a tuple of Ligand objects

Arguments:
	filename::str
		Path to the file holding the ligand data
	
Returns:
	ligands::tuple(Ligand)
		A tuple of Ligand objects of length n, n being the number of ligands
 ```
 
  
 
tanimoto_distance(bit_vec1, bit_vec2):
```
Calculates the Tanimoto distance for two bit vectors of the same length. Tanimoto coefficient is defined as intersection(a,b)/union(a,b),
and Tanimoto distance is defined as 1 - Tanimoto Coefficient

Arguments:
	bit_vec1::array(bool/int)
		A numpy array of length n, in the case of these molecular fingerprint scores n=1024. Each entry is zero or one based on the absence or presence of 		   the motif that each feature represenets. 
	bit_vec2::array(bool/int)
		A numpy array of length n, in the case of these molecular fingerprint scores n=1024. Each entry is zero or one based on the absence or presence of 	           the motif that each feature represenets.
	
Returns:
	distance::float
		Tanimoto distance as defined in the function description between bit_vec1 and bit_vec2
```
