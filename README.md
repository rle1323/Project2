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

clustering_similarity(labels1, labels2)
```
Calculates the Jaccard similarity between two sets of clustering results. Inspired by the R function linked here:
https://rdrr.io/cran/clusteval/man/jaccard_indep.html

Arguments:
	labels1::[int]
		Prediction of cluster memberships from first clustering
	labels2::[int]
		Prediction of cluster memberships from second memberships
	
Returns:
	similarity::float
		Jaccard similarity defined as n_11/(n_11 + n_10 + n_01), where n_11 is number of observation pairs that are comembers in both clusterings, 
		n_10 is number of observation pairs that are comembers in first cluster but not the second, and n_01 is number of observatoin pairs that are 
		comembers in second cluster but not the first
```

build_distance_matrix(ligands):
```
Builds an nxn distance matrix for a set of n ligands, where distance is defined as (1-Tanimoto coefficient)

Arguments:
	ligands::[Ligand]
		List of n Ligand objects that the distance matrix will be made for. Function unpacks the relevant attributes from each Ligand in the list
	
Returns:
	distance_matrix::array(float)
		An nxn distance matrix where each entry distance_matrix[i,j] corresponds to the Tanimoto distance between Ligand i and Ligand j
		in the input list
```
silhouette_score(ligands, labels, distance_matrix=None):
```
Calculates the mean silhouette score from the results of a clustering. This score can be used as a general quality metric of a clustering.
Silhouette score is defined as (b-a)/max(a,b), where a is the mean distance from a point to the other points in its clustser,
and b is defined as the mean distance from a point to the points in the next closest cluster. Has range[-1,1]. 

Arguments:
	ligands::[Ligand]
		List of n ligand objects that have been clustered
	labels::[int]
		List of n cluster labels for each of the n ligands
	(Optional) distance_matrix::np.array(int)
		Distance matrix of the supplied ligands. If not provided, the function calls build_distance_matrix() to build the matrix
	
Returns:
	silhouette_score::float
		Mean silhouette score of all of the ligands in this clustering.
```

### Classes

Ligand:
```
A class to hold and transform data provided in ligand_information.csv
```
__init__(self, id, score, smiles, on_bits):
```
Provided a ligand ID, Vina score, SMILES string, and dense molecular fingerprint vector, initializes a Ligand object and transforms the dense 
molecular fingerprint into a bit vector.

Arguments:
	id::int
		Unique ID of the ligand in the dataset
	score::float
		Vina score of the ligand against ACE-Spike2 interface
	smiles::string
		SMIlES string of the molecule, holding structural and composition encoding
	on_bits::string
		Extended Connectivity Fingerprint (ECFP) of the compound
	
Returns:
	None
```
fingerprint_array(self):
```
Transforms the ECFP into a bit vector and stores the bit vector in a new attribute

Arguments:
	None
		
Returns:
	None
```

Cluster:
```
A simple class holding information for the clusters in either Kmeans or agglomerative clustering
```
__init__(self, members, label, centroid=None):
```
Takes in a list of ligands as members, a label for the cluster, and an optional centroid argument (for Kmeans) and initializes a Cluster object
```
