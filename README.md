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

__Ligand__:
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

__Cluster__:
```
A simple class holding information for the clusters in either Kmeans or agglomerative clustering
```
init(self, members, label, centroid=None):
```
Takes in a list of ligands as members, a label for the cluster, and an optional centroid argument (for Kmeans) and initializes a Cluster object
```

__Clustering__:
```
A parent class for each of the two clustering methods, holding shared attributes and methods
```

init(self, num_clusters, seed=2000):
```
Initializes a Clustering object and optionally sets a seed.

Arguments:
	num_clusters::int
		Number of desired clusters that the data will be placed into
	seed::int
		Random seed for reproducibility
		
Returns:
	None
```
init_clusters(self):
```
Generic method to set initial locations of the clusters in the feature space. Specific methods are in the child class definitions. 
```
cluster(self):
```
Generic method that performs the specific clustering method implemented in the child classes
```

__HierarchicalClustering__:
```
Complete implementation of agglomerative hierarchical clustering. To use this class, initialize a HierarchicalClustering object and use it to call the 
cluster() method on a set of ligands. Child class of Clustering
```
init(self, num_clusters, seed):
```
Uses the same attributes as the parent class Clustering, and initializes a HierarchicalClustering object
Arguments:
	num_clusters::int
		Number of desired clusters that the data will be placed into
	(Optional) seed::int
		Random seed for reproducibility

Returns:
	None
```
build_id_dictionary(self, ligands):
```
Builds a dictionary with ligand ids as keys and the ligand row (0-indexed) in the dataset as the value. Saves the dictionary as an attribute

Arguments:
	ligands::[Ligand]
		List of ligand objects for which the dictionary will be built
		
Returns:
	None
```

init_clusters(self, ligands):
```
An agglomerative hierarchical method is used for this clustering, so every ligand is first assigned to its own cluster. A list of Cluster objects
is set as a class attribute

Arguments:
	ligands::[Ligand]
		List of ligand objects on which the clustering is being performed
Returns:
	None
```
complete_linkage(self, cluster1, cluster2): 
```
Finds the complete linkage betweeen two clusters. Complete linkage is defined as the greatest distance between any two members of the two clusters.

Arguments:
	cluster1::Cluster
		First cluster whose complete linkage with cluster2 is being assessed
	cluster2::Cluster
		Second cluster whose complete linkage with cluster1 is being assessed
		
Returns:
	complete_linkage::float
		Complete linkage between cluster1 and cluster2
```
nearest_neighboring_clusters(self):
```
Finds the two closest clusters as defined by the linkage criterion, returns their indices in self.clusters

Arguments:
	None
		
Returns: 
	nearest_neighbor1::int
		Index of the first of the two nearest clusters. This index marks where this neighbor is held in the classes self.clusters list
	nearest_neighbor2::int
		Index of the second of the two nearest clusters. This index marks where this neighbor is held in the classes self.clusters list.
```

merge_nearest_clusters(self, i, j):
```
Merges two clusters, found at positions i and j in the self.clusters list. Performs merging by iteratively adding the members of 
self.clusters[j] to the member list of self.clusters[i] and removing self.clusters[j] from the list of clusters

Arguments:
	i::int
		Position of the cluster that is growing in membership in self.clusters
	j::int
		Position of the cluster whose members are merged into another cluster and eventually removed
		
Returns:
	None
```
cluster(self, ligands, distance_matrix=None):
```
Wrapper method for hierarchical clustering that performs the above methods in their proper order.

Arguments:
	ligands::[Ligand]
		List of Ligand objects that are being clustered
	(Optional) distance_matrix::array(float)
		A distance matrix of the ligands passed into this method. This argument is optional, and this method will build the 
		distance matrix if it is not specified
		
Returns:
	labels::[int]
		List of labels with the same length as the input ligand list. Each label corresponds to the ligand at that same index in ligands
```
__PartitionClustering__:
```
Complete implementation of Kmeans partition clustering. To use this class, initialize a PartitionClustering object and use it to call the 
cluster() method on a set of ligands. Child class of Clustering.
```

init(self, num_clusters, seed=1998, max_iterations=1000):
```
Uses the same attributes as the parent class Clustering with an additional "max_iterations" attribute, and initializes a PartitionClustering object

Arguments:
	num_clusters::int
		Number of desired clusters that the data will be placed into
	(Optional) seed::int
		Random seed for reproducibility
	(Optional) max_iterations::int
		Number of iterations to perform if convergence is not reached
		
Returns:
	None
```
init_clusters(self, ligands):
```
Uses a simplified kmeans++ initialization scheme to set initial clusters with centroid locations

Arguments:
	ligands::[Ligand]
		List of ligand objects on which the clustering is being performed
		
Returns:
	None
```
assign_cluster_membership(self, ligands):
```
Assigns each ligand to the nearest cluster as defined by Tanimoto distance between the ligand and cluster centroids.
Adds the ligand to the cluster's memmbership list

Arguments:
	ligands::[Ligand]
		List of ligand objects on which the clustering is being performed
		
Returns:
	None
```
 update_cluster_centroids(self):
```
Recomputes the centroid after cluster members are changed. Each centroid feature is the mode of the feature in the cluster's members

Arguments:
	None
		
Returns:
	None
```
check_convergence(self, old_memberships):
```
Checks the convergence condition for partition clustering. That is, checks whether cluster membership has changed for any of the ligands
from one iteration to the next

Arguments:
	old_memberships::[[Ligand]]
		List of membership lists, each of which contains the ligands belonging to each cluster
		
Returns:
	::bool
		True if the convergence criterion is met, false otherwise
```
cluster(self, ligands):
```
Wrapper function for partition clustering. Partition clustering is implemented with Kmeans++ algorithm

Arguments:
	ligands::[Ligand] 
		The list of ligands that are being clustered 

Returns:
	labels::[int]
		List of labels with the same length as the input ligand list. Each label corresponds to the ligand at that same index in ligands
```
