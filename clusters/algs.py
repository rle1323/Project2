import numpy as np 
import pandas as pd
import copy as cp
import statistics as stats

############### Function Definitions ###############

def read_ligand_data(filename):
	"""
	Reads in the provided ligand data from csv format, and stores the ligands as a tuple of Ligand objects

	Arguments:
		filename::str
			Path to the file holding the ligand data
	
	Returns:
		ligands::tuple(Ligand)
			A tuple of Ligand objects of length n, n being the number of ligands
	"""

	ligands = []
	# read the csv file storing ligand information into a panda frame
	file_pd = pd.read_csv(filename)

	# for each entry in the file, create a Ligand object and put that object into a list of the ligands
	for i, row in file_pd.iterrows():
		new_ligand = Ligand(
			id=row.loc["LigandID"], 
			score=row.loc["Score"], 
			smiles=row.loc["SMILES"], 
			on_bits=row.loc["OnBits"]
		)
		ligands.append(new_ligand)
	
	ligands = tuple(ligands)

	return ligands

def tanimoto_distance(bit_vec1, bit_vec2):
	"""
	Calculates the Tanimoto distance for two bit vectors of the same length. Tanimoto coefficient is defined as intersection(a,b)/union(a,b),
	and Tanimoto distance is defined as 1 - Tanimoto Coefficient

	Arguments:
		bit_vec1::array(bool/int)
			A numpy array of length n, in the case of these molecular fingerprint scores n=1024. Each entry is zero or one based on the absence or presence of the 
			motif that each feature represenets. 
		bit_vec2::array(bool/int)
			A numpy array of length n, in the case of these molecular fingerprint scores n=1024. Each entry is zero or one based on the absence or presence of the 
			motif that each feature represenets.
	
	Returns:
		distance::float
			Tanimoto distance as defined in the function description between bit_vec1 and bit_vec2.
	"""

	# set the bit vecs as sets so intersection and union can be calculated
	temp1 = set(np.flatnonzero(bit_vec1))
	temp2 = set(np.flatnonzero(bit_vec2))
	
	# find intersection and union
	intersection = temp1 & temp2
	union = temp1 | temp2
	
	similarity = len(intersection)/len(union)
	distance = 1 - similarity
		
	return distance


def clustering_similarity(labels1, labels2):
	"""
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
	"""
	n_11 = 0
	n_other = 0

	for i in range(len(labels1)):
		for j in range(len(labels1))[i+1:]:
			# check equality pairs in both sets of results
			result1 = (labels1[i] == labels1[j])
			result2 = (labels2[i] == labels2[j])
			
			if result1 == result2 and result1 == True:
				n_11 += 1
			elif result1 != result2:
				n_other += 1
	
	similarity = n_11/(n_11 + n_other)

	return similarity


def build_distance_matrix(ligands):
	"""
	Builds an nxn distance matrix for a set of n ligands, where distance is defined as (1-Tanimoto coefficient)

	Arguments:
		ligands::[Ligand]
			List of n Ligand objects that the distance matrix will be made for. Function unpacks the relevant attributes from each Ligand in the list
	
	Returns:
		distance_matrix::array(float)
			An nxn distance matrix where each entry distance_matrix[i,j] corresponds to the Tanimoto distance between Ligand i and Ligand j
			in the input list
	"""
	num_ligands = len(ligands)
	distance_matrix = np.zeros([num_ligands, num_ligands])

	# fill out distance to each other ligand for each ligand. set a ligand's self-distance to infinite so that it doesnt come up as "closest neighbor"
	for i,ligand1 in enumerate(ligands): 
		for j,ligand2 in enumerate(ligands[i+1:]):
			jdx = j+i+1
			distance_matrix[i,jdx] = tanimoto_distance(ligand1.bit_vector, ligand2.bit_vector)
		
	# fill in bottom of matrix, and set diagonal to infinite
	distance_matrix = distance_matrix + distance_matrix.T
	np.fill_diagonal(distance_matrix, float("inf"))

	return distance_matrix

def silhouette_score(ligands, labels, distance_matrix=None):
	"""
	Calculates the mean silhouette score from the results of a clustering. This score can be used as a general quality metric of a clustering.
	Silhouette score is defined as (b-a)/max(a,b), where a is the mean distance from a point to the other points in its clustser,
	and b is defined as the mean distance from a point to the points in the next closest cluster. Has range[-1,1]. 

	Arguments:
		ligands::[Ligand]
			List of n ligand objects that have been clustered
		labels::[int]
			List of n cluster labels for each of the n ligands
		distance_matrix::np.array(int) or None
			Distance matrix of the supplied ligands. If not provided, the function calls build_distance_matrix() to build the matrix
	
	Returns:
		silhouette_score::float
			Mean silhouette score of all of the ligands in this clustering.
	"""
	# if a distance matrix is not supplied, compute it
	if distance_matrix is None:
		distance_matrix = build_distance_matrix(ligands)
		
	num_clusters = max(labels)
	
	# generate a score for every datapoint
	point_scores = []
	for i in range(len(labels)):
		ligand1 = ligands[i]
		label1 = labels[i]
		cluster_distances = [[] for _ in range(num_clusters)]

		for j in range(len(labels)):
			ligand2 = ligands[j]
			label2 = labels[j]

			# don't compare a ligand to itself for this calculation
			if ligand1 is not ligand2:
				# get distance between the two ligands, and put the distance in the appropriate sub-list
				dist = distance_matrix[i,j]
				cluster_distances[label2-1].append(dist)
			# after distances have been indexed between this ligand and every other, calculate the means of the sub-lists

		# if a cluster is empty, fill its distance list with dummy value
		for i in range(num_clusters):
			if not cluster_distances[i]:
				cluster_distances[i] = [2]
		
		# calculate means of distance from current ligand to each cluster's members, including its own cluster
		means = [stats.mean(sub) for sub in cluster_distances]
		intra_cluster_mean = means[label1-1]
		
		# find the mean distance of the next nearest cluster
		nearest_inter_cluster_mean = 2
		for x in range(num_clusters):
			if x == (label1-1): continue
			else:
				if means[x] < nearest_inter_cluster_mean:
					nearest_inter_cluster_mean = means[x] 
			
		# do final calculation of silhouette score, which equals (b-a)/max(a,b)
		silhouette = (nearest_inter_cluster_mean - intra_cluster_mean)/max(nearest_inter_cluster_mean, intra_cluster_mean)
		point_scores.append(silhouette)
		
	silhouette_score = stats.mean(point_scores)
		
	return silhouette_score

############### Class Definitions ###############

class Ligand():
	"""
	A class to hold and transform data provided in ligand_information.csv
	"""

	def __init__(self, id, score, smiles, on_bits):
		"""
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
		"""
		self.id = id
		self.score = score
		self.smiles = smiles
		self.on_bits_string = on_bits
		self._fingerprint_array()


	def _fingerprint_array(self):
		"""
		Transforms the ECFP into a bit vector and stores the bit vector in a new attribute

		Arguments:
			None
		
		Returns:
			None
		"""

		# split the string into separate elements, and convert them to integers
		dense_vec = [int(i) for i in self.on_bits_string.split(",")]

		# fill the entire vector with zeros, then put 1s at the select indices 
		self.bit_vector = np.zeros(1024)
		self.bit_vector[dense_vec] = 1


class Cluster():
	"""
	A simple class holding information for the clusters in either Kmeans or agglomerative clustering
	"""

	def __init__(self, members, label, centroid=None):
		"""
		Takes in a list of ligands as members, a label for the cluster, and an optional centroid argument (for Kmeans) and initializes a Cluster object
		"""
		self.members = members
		self.label = label
		self.centroid = centroid


class Clustering():
	"""
	A parent class for each of the two clustering methods, holding shared attributes and methods
	"""

	def __init__(self, num_clusters, seed=2000):
		"""
		Initializes a Clustering object and optionally sets a seed.

		Arguments:
			num_clusters::int
				Number of desired clusters that the data will be placed into
			seed::int
				Random seed for reproducibility
		
		Returns:
			None
		"""
		np.random.seed(seed)
		self.clusters = []
		self.num_clusters = num_clusters
	
	def _init_clusters(self):
		"""
		Generic method to set initial locations of the clusters in the feature space. Specific methods are in the child class definitions. 
		"""
		pass
	
	def cluster(self):
		"""
		Generic method that performs the specific clustering method implemented in the child classes
		"""
		pass
		


class HierarchicalClustering(Clustering):
	"""
	Complete implementation of agglomerative hierarchical clustering. To use this class, initialize a HierarchicalClustering object and use it to call the 
	cluster() method on a set of ligands. Child class of Clustering
	"""

	def __init__(self, num_clusters, seed=1998):
		"""
		Uses the same attributes as the parent class Clustering, and initializes a HierarchicalClustering object
		Arguments:
			num_clusters::int
				Number of desired clusters that the data will be placed into
			(Optional) seed::int
				Random seed for reproducibility
		
		Returns:
			None
		"""
		super().__init__(num_clusters, seed)
	
	def _build_id_dictionary(self, ligands):
		"""
		Builds a dictionary with ligand ids as keys and the ligand row (0-indexed) in the dataset as the value. Saves the dictionary as an attribute

		Arguments:
			ligands::[Ligand]
				List of ligand objects for which the dictionary will be built
		
		Returns:
			None
		"""
		self.ligand_dictionary = {}
		for i, ligand in enumerate(ligands):
			self.ligand_dictionary[ligand.id] = i

	
	def _init_clusters(self, ligands):
		"""
		An agglomerative hierarchical method is used for this clustering, so every ligand is first assigned to its own cluster. A list of Cluster objects
		is set as a class attribute

		Arguments:
			ligands::[Ligand]
				List of ligand objects on which the clustering is being performed
		"""
		for i, ligand in enumerate(ligands):
			new_cluster = Cluster(members=[ligand], label = i+1)
			self.clusters.append(new_cluster)

	
	def _complete_linkage(self, cluster1, cluster2): 
		"""
		Finds the complete linkage betweeen two clusters. Complete linkage is defined as the greatest distance between any two members of the two clusters.

		Arguments:
			cluster1::Cluster
				First cluster whose complete linkage with cluster2 is being assessed
			cluster2::Cluster
				Second cluster whose complete linkage with cluster1 is being assessed
		
		Returns
			complete_linkage::float
				Complete linkage between cluster1 and cluster2
		"""
		distances = []
		for ligand1 in cluster1.members:
			i = self.ligand_dictionary[ligand1.id]
			for ligand2 in cluster2.members:
				j = self.ligand_dictionary[ligand2.id]
				distances.append(self.distance_matrix[i,j])
		
		# get maximum distance between any two members
		complete_linkage = max(distances)
		return complete_linkage


	def _nearest_neighboring_clusters(self):
		"""
		Finds the two closest clusters as defined by the linkage criterion, returns their indices in self.clusters

		Arguments:
			None
		
		Returns: 
			nearest_neighbor1::int
				Index of the first of the two nearest clusters. This index marks where this neighbor is held in the classes self.clusters list
			nearest_neighbor2::int
				Index of the second of the two nearest clusters. This index marks where this neighbor is held in the classes self.clusters list.
		"""
		# initial a matrix of linkage distances
		linkage_distances = np.ones([len(self.clusters), len(self.clusters)])
		np.fill_diagonal(linkage_distances, float("inf"))

		# find each cluster's linkage with every other cluster
		for i, cluster1 in enumerate(self.clusters):
			for j, cluster2 in enumerate(self.clusters[i+1:]):
				jdx = i+j+1
				linkage_distance = self._complete_linkage(cluster1, cluster2)
				if linkage_distance == 0:
					return i,jdx
				else:
					linkage_distances[i,jdx] = linkage_distance
		
		# find the minimum distance among the clusters
		nearest_neighbor1, nearest_neighbor2 = np.unravel_index(linkage_distances.argmin(), linkage_distances.shape)
		return nearest_neighbor1, nearest_neighbor2
		

	def _merge_nearest_clusters(self, i, j):
		"""
		Merges two clusters, found at positions i and j in the self.clusters list. Performs merging by iteratively adding the members of 
		self.clusters[j] to the member list of self.clusters[i] and removing self.clusters[j] from the list of clusters

		Arguments:
			i::int
				Position of the cluster that is growing in membership in self.clusters
			j::int
				Position of the cluster whose members are merged into another cluster and eventually removed
		
		Returns:
			None
		"""
		# add the members of the second cluster to the membership list of the first cluster
		for member in self.clusters[j].members:
			self.clusters[i].members.append(member)

		# remove the 2nd cluster from the cluster list
		del self.clusters[j]


	
	def cluster(self, ligands, distance_matrix=None):
		"""
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
		"""
		# make id dictionary for the ligands
		self._build_id_dictionary(ligands)

		# if a distance matrix is not supplied, make it
		if distance_matrix is None:
			self.distance_matrix = build_distance_matrix(ligands)
		else:
			self.distance_matrix = distance_matrix

		# initialize clusters
		self._init_clusters(ligands)
		cluster_count = len(self.clusters)

		#continue until we have the desired number of clusters
		while cluster_count > self.num_clusters:
			
			# get indices of two nearest clusters according to linkage criterion
			nearest1, nearest2 = self._nearest_neighboring_clusters()
			
			# merge these two clusters
			self._merge_nearest_clusters(nearest1, nearest2)
			cluster_count -= 1

		# fix cluster labels so that they are of range(1, num_clusters)
		for i in range(1, len(self.clusters)+1):
			self.clusters[i-1].label = i

		# go through each ligand, find which cluster it is in, and save that label
		labels = []
		for ligand in ligands:
			for cluster in self.clusters:
				if ligand in cluster.members:
					labels.append(cluster.label)

		return labels




class PartitionClustering(Clustering):
	"""
	Complete implementation of Kmeans partition clustering. To use this class, initialize a PartitionClustering object and use it to call the 
	cluster() method on a set of ligands. Child class of Clustering.
	"""

	def __init__(self, num_clusters, seed=1998, max_iterations=1000):
		"""
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
		"""
		super().__init__(num_clusters, seed)
		self.max_iterations = max_iterations
	

	def _init_clusters(self, ligands):
		"""
		Uses a simplified kmeans++ initialization scheme to set initial clusters with centroid locations

		Arguments:
			ligands::[Ligand]
				List of ligand objects on which the clustering is being performed
		
		Returns:
			None
		"""
		# sample a random datapoint to serve as first cluster centroid, make a new cluster object and add it to the list
		centroid1 = np.random.choice(ligands).bit_vector
		cluster1 = Cluster(members=[], label=1, centroid=centroid1)
		self.clusters.append(cluster1)

		# find where the remaining initial centroids will be
		for i in range(self.num_clusters - 1):
			
			# list to store distances of point from nearest centroid
			distances = []
			for ligand in ligands:
				# set an initial distance value that is higher than maximum tanimoto distance (1)
				min_dist = 2

				# compute distance from current ligand to all previously set centroids
				for cluster in self.clusters:
					calculated_dist = tanimoto_distance(ligand.bit_vector, cluster.centroid)
					min_dist = min(calculated_dist, min_dist)
				distances.append(min_dist)
			
			# select data point with maximum distance from it's closest centroid neighbor as next centroid
			new_centroid_ligand = ligands[np.argmax(np.array(distances))]
			new_cluster = Cluster(members=[], label=i+2, centroid=new_centroid_ligand.bit_vector)
			self.clusters.append(new_cluster)


	def _assign_cluster_membership(self, ligands):
		"""
		Assigns each ligand to the nearest cluster as defined by Tanimoto distance between the ligand and cluster centroids.
		Adds the ligand to the cluster's memmbership list

		Arguments:
			ligands::[Ligand]
				List of ligand objects on which the clustering is being performed
		
		Returns:
			None
		"""
		# whipe the members of each cluster
		for i in range(len(self.clusters)):
			self.clusters[i].members = []

		# for each ligand, find closest cluster centroid and add the ligand to that cluster's member list
		for ligand in ligands:
			distances = []
			for i, cluster in enumerate(self.clusters):
				dist = tanimoto_distance(cluster.centroid, ligand.bit_vector)
				distances.append(dist)

			nearest = np.argmin(np.array(distances))
			# if this is the first member of the cluster, make a new list for the nearest cluster. otherwise, append the ligand to the member list of the nearest cluster
			if not self.clusters[nearest].members:
				self.clusters[nearest].members = [ligand]
			else:
				self.clusters[nearest].members.append(ligand)
			

	def _update_cluster_centroids(self):
		"""
		Recomputes the centroid after cluster members are changed. Each centroid feature is the mode of the feature in the cluster's members

		Arguments:
			None
		
		Returns:
			None
		"""
		bit_space_size = self.clusters[0].members[0].bit_vector.shape[0]

		for i, cluster in enumerate(self.clusters):
			# for each feature, find the mode value of the cluster's members and set that to the centroid feature
			member_bit_vecs = np.array([ligand.bit_vector for ligand in cluster.members])

			# make a new centroid array
			new_centroid = np.zeros(bit_space_size)
			if len(cluster.members) == 0:
				continue
			elif len(cluster.members) == 1:
				new_centroid = cluster.members[0]		
			else:
				for j in range(bit_space_size):
					feature = member_bit_vecs[:,j]
					mode = round(np.mean(feature))
					new_centroid[j] = mode
			
			self.clusters[i].centroid = new_centroid
	

	def _check_convergence(self, old_memberships):
		"""
		Checks the convergence condition for partition clustering. That is, checks whether cluster membership has changed for any of the ligands
		from one iteration to the next

		Arguments:
			old_memberships::[[Ligand]]
				List of membership lists, each of which contains the ligands belonging to each cluster
		
		Returns:
			::bool
				True if the convergence criterion is met, false otherwise
		"""
		for i, cluster in enumerate(self.clusters):
			old_membership = old_memberships[i]
			old_ids = np.array([ligand.id for ligand in old_membership])
			
			# retrieve the ids of the newly assigned members
			new_ids = np.array([ligand.id for ligand in cluster.members])
			# check if old memberships are equal to the new ones
			if not np.array_equal(old_ids, new_ids):
				return False
		
		return True
		

	def cluster(self, ligands):
		"""
		Wrapper function for partition clustering. Partition clustering is implemented with Kmeans++ algorithm

		Arguments:
			ligands::[Ligand] 
				The list of ligands that are being clustered 

		Returns:
			labels::[int]
				List of labels with the same length as the input ligand list. Each label corresponds to the ligand at that same index in ligands
		"""
		# perform cluster initialization with kmeans++
		self._init_clusters(ligands)
		
		# perform first cluster assignments
		self._assign_cluster_membership(ligands)

		for i in range(self.max_iterations): 
			if i > 0:
				# save the old memberships so that we can check for convergence later
				old_memberships = [cp.deepcopy(cluster.members) for cluster in self.clusters]

				# reassign ligands to nearest cluster
				self._assign_cluster_membership(ligands)

				# check for convergence, and return if converged
				if self._check_convergence(old_memberships) == True: 
					print("Converged after ", i, "iterations")
					break

				# update the cluster centroids
				self._update_cluster_centroids()
				
			else: # if this is the first iteration, just set the cluster centroids and move on
				self._update_cluster_centroids()
		
		# go through each ligand, find which cluster it is in, and save that label
		labels = []
		for ligand in ligands:
			for cluster in self.clusters:
				if ligand in cluster.members:
					labels.append(cluster.label)

		return labels