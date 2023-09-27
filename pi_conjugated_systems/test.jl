using ClusterMeanField
using RDM
using QCBase
using InCoreIntegrals
using PyCall
using ActiveSpaceSolvers
using LinearAlgebra
using Printf
using NPZ
using JLD2
clusters = [[1,2,6,7,8],[3,4,5]]
init_fspace = [(3,3),(1,1)]


rdm1 = RDM1(n_orb(ints))

# define clusters
# clusters = [MOCluster(i,collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)



norb=8
ansatze = [FCIAnsatz(5, 3, 3),FCIAnsatz(3,1,1)]
proj_vec = ClusterMeanField.projection_vector(ansatze, norb)
display(proj_vec)
projection_vec=ClusterMeanField.projection_vector2(ansatze, norb,clusters)
println(projection_vec)
display(projection_vec)