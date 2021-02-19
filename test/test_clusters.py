import pytest
import numpy as np
from clusters import algs

# read in the ligand dataset because it will be needed
ligands = algs.read_ligand_data("ligand_information.csv")

# tests that the readin returns a list of the proper length
def test_readin():
	assert(len(ligands) == 8524)


# tests my tanimoto distance function
def test_tanimoto():
	# these two should have distance of 1
	bit_vec1=np.array([1,1,0,0])
	bit_vec2=np.array([0,0,1,1])
	assert(algs.tanimoto_distance(bit_vec1, bit_vec2) == 1)

	# these two should have distance of zero
	bit_vec1=np.array([1,1,0,0])
	bit_vec2=np.array([1,1,0,0])
	assert(algs.tanimoto_distance(bit_vec1, bit_vec2) == 0)

	# these two should have distance of .5
	bit_vec1=np.array([1,0,0,0])
	bit_vec2=np.array([1,0,1,0])
	assert(algs.tanimoto_distance(bit_vec1, bit_vec2) == .5)


def test_clustering_similarity():
	# testing the output of my function with some dummy labels directly against an R implementation of my function described here: https://rdrr.io/cran/clusteval/man/jaccard_indep.html
	# these should have similarity 1
	labels1 = [1,1,2,2]
	labels2 = [1,1,2,2]
	assert(algs.clustering_similarity(labels1,labels2) == 1)

	# these should have similarity .5
	labels1 = [1,1,1,1]
	labels2 = [1,1,1,2]
	assert(algs.clustering_similarity(labels1,labels2) == .5)

	# these should have similarity .2
	labels1 = [1,1,2,1]
	labels2 = [1,2,1,1]
	assert(algs.clustering_similarity(labels1,labels2) == .2)


def test_build_distance_matrix():
	# makes sure that the dimensionality of the distance matrix is correct
	distance_mat = algs.build_distance_matrix(ligands[:500])
	assert(distance_mat.shape == (500,500))


def test_silhouette_score():
	# tests that the silhouette score falls between -1 and 1
	fake_labels = [1 for i in range(250)]
	for i in range(250):
		fake_labels.append(2)
	fake_dist_mat = algs.build_distance_matrix(ligands[:500])
	score = algs.silhouette_score(ligands[:500], fake_labels, fake_dist_mat)
	assert(score >= -1 and score <= 1)

def test_partitioning():
	# gets the clustering results from Kmeans clustering, and checks that the number of different labels is correct
	test_cluster5 = algs.PartitionClustering(5, seed = 2)
	labels5 = test_cluster5.cluster(ligands[:500])
	assert(1 in labels5 and 2 in labels5 and 3 in labels5 and 4 in labels5 and 5 in labels5)

	test_cluster2 = algs.PartitionClustering(10, seed = 6)
	labels2 = test_cluster2.cluster(ligands[:500])
	assert(1 in labels2 and 2 in labels2)

	# for each ligand, check that the closest cluster centroid is the cluster that it belongs to
	test_cluster2 = algs.PartitionClustering(10, seed = 6)
	labels2 = test_cluster2.cluster(ligands[:500])
	for i, ligand in enumerate(ligands[:500]):
		distances = []
		for cluster in test_cluster2.clusters:
			distances.append(algs.tanimoto_distance(ligand.bit_vector, cluster.centroid))
		assert(np.argmin(np.array(distances)) + 1 == labels2[i])


def test_hierarchical():
	# gets the clustering results from agglomerative clustering, and checks that the number of different labels is correct
	# use tiny datasets because this ish takes forever
	test_cluster5 = algs.HierarchicalClustering(5)
	labels5 = test_cluster5.cluster(ligands[:100])
	assert(1 in labels5 and 2 in labels5 and 3 in labels5 and 4 in labels5 and 5 in labels5)

	test_cluster2 = algs.HierarchicalClustering(10, seed = 6)
	labels2 = test_cluster2.cluster(ligands[:100])
	assert(1 in labels2 and 2 in labels2)

