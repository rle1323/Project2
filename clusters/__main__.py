import numpy as np 
from clusters import algs
import umap
import matplotlib.pyplot as plt


# read in the ligand data 
ligands = algs.read_ligand_data("./ligand_information.csv")


################## QUESTION 2 ##################
'''
# get just the bit vectors of the ligands and save it as a list of lists
ligand_bit_vecs = []
for ligand in ligands:
    bv = []
    for bit in ligand.bit_vector:
        bv.append(bit)
    ligand_bit_vecs.append(bv)
# generate a umap object and reduce the data
umap_reduction = umap.UMAP()
umap_embedding = umap_reduction.fit_transform(ligand_bit_vecs)

# plot the embedded ligand data in 2d
plt.scatter(umap_embedding[:, 0], umap_embedding[:, 1], 5)
plt.show()
'''

################## QUESTION 3 ##################
'''
# precompute the distance matrix so that it doesnt need to be recomputed every time
print("Building distance mat")
#distance_matrix = algs.build_distance_matrix(ligands)
print("distance mat built")
# set up a kmeans clustering scheme to determine which value of k performs the best. try values of k between 2-15
k_vals = [i for i in range(2, 16)]
scores = []
for k in k_vals:
    test_cluster = algs.PartitionClustering(num_clusters = k, seed = 6, max_iterations=25)
    labels = test_cluster.cluster(ligands)
    score = algs.silhouette_score(ligands, labels, distance_matrix)
    print(score)
    scores.append(score)
# make the line plot of k vs silhouette score
plt.plot(k_vals, scores)
plt.title('Silhouette score by K -- Partition Clustering')
plt.xlabel('K')
plt.ylabel('Score')
plt.show()

# get just the bit vectors of the ligands and save it as a list of lists
ligand_bit_vecs = []
for ligand in ligands:
    bv = []
    for bit in ligand.bit_vector:
        bv.append(bit)
    ligand_bit_vecs.append(bv)
# generate a umap object and reduce the data
umap_reduction = umap.UMAP(transform_seed = 3)
umap_embedding = umap_reduction.fit_transform(ligand_bit_vecs)

# generate a clustering with k = 2
test_cluster = algs.PartitionClustering(num_clusters = 2, seed = 6, max_iterations=25)
labels = test_cluster.cluster(ligands)

plt.scatter(umap_embedding[:, 0], umap_embedding[:, 1], 5, c = labels)
plt.show()



################## QUESTION 5 ##################

# generate random subset of 2000 ligands
np.random.seed(5)
num_ligands = (len(ligands))
random_indices = np.random.choice(num_ligands, size=2000, replace=False)
ligand_subset = [ligands[idx] for idx in random_indices]

# get just the bit vectors of the ligands and save it as a list of lists
ligand_bit_vecs = []
for ligand in ligand_subset:
    bv = []
    for bit in ligand.bit_vector:
        bv.append(bit)
    ligand_bit_vecs.append(bv)
# generate a umap object and reduce the data
umap_reduction = umap.UMAP(transform_seed = 3)
umap_embedding = umap_reduction.fit_transform(ligand_bit_vecs)


# generate a clustering
test_cluster = algs.HierarchicalClustering(num_clusters=5)
labels = test_cluster.cluster(ligand_subset)

# show plot
plt.scatter(umap_embedding[:, 0], umap_embedding[:, 1], 5, c = labels)
plt.show()


################## QUESTION 7+8 ##################

# generate random subset of 2000 ligands
np.random.seed(5)
num_ligands = (len(ligands))
random_indices = np.random.choice(num_ligands, size=2000, replace=False)
ligand_subset = [ligands[idx] for idx in random_indices]

# generate distance matrix for subset
subset_distance_matrix = algs.build_distance_matrix(ligand_subset)

# get k-means clustering
partition_cluster = algs.PartitionClustering(num_clusters = 5, seed = 6, max_iterations=25)
partition_labels = partition_cluster.cluster(ligand_subset)
print("PARTITION CLUSTERING SCORE", algs.silhouette_score(ligand_subset, partition_labels, subset_distance_matrix))

# get hierarchical clustering
hierarchical_cluster = algs.HierarchicalClustering(num_clusters=5)
hierarchical_labels = hierarchical_cluster.cluster(ligand_subset,  distance_matrix = subset_distance_matrix)

# get quality scores for each
print("AGGLOM CLUSTERING SCORE", algs.silhouette_score(ligand_subset, hierarchical_labels, subset_distance_matrix))

# get clustering similarity 
print("CLUSTERING SIMILARITY", algs.clustering_similarity(partition_labels, hierarchical_labels) )
'''


################## QUESTION 9 ##################
'''
# get k-means clustering
partition_cluster = algs.PartitionClustering(num_clusters = 5, seed = 6, max_iterations=25)
partition_labels = partition_cluster.cluster(ligands)

# group the scores for each cluster
cluster_scores = [[], [], [], [], []]

for ligand, label in zip(ligands, partition_labels):
    cluster_scores[label-1].append(ligand.score)

for i, scores in enumerate(cluster_scores):
    plt.hist(scores)
    plt.xlabel("Vina Score")
    plt.title("Distribution of Vina Scores for cluster "+ str(i+1))
    plt.show()
'''

################## QUESTION 10 ##################
# get k-means clustering
partition_cluster = algs.PartitionClustering(num_clusters = 5, seed = 6, max_iterations=25)
partition_labels = partition_cluster.cluster(ligands)

# find the 
# group the scores for each cluster
clusters = [[], [], [], [], []]

for ligand, label in zip(ligands, partition_labels):
    clusters[label-1].append(ligand)

# find the top scoring ligand in each cluster and print
for i, cluster in enumerate(clusters):
    max_score = float("-inf")
    for ligand in cluster:
        if ligand.score > max_score:
            max_score = ligand.score
            highest_scoring = ligand
    print("Highest scoring in cluster", i+1, "is ligand", ligand.id, "with Vina score", ligand.score)


