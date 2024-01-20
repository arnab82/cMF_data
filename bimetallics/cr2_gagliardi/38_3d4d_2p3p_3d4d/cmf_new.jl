using QCBase
using ClusterMeanField 
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using ActiveSpaceSolvers

h0 = npzread("/Users/ayush/workspace/cMF_data/bimetallics/cr2_gagliardi/38_3d4d_2p3p_3d4d/ints_h0.npy")
h1 = npzread("/Users/ayush/workspace/cMF_data/bimetallics/cr2_gagliardi/38_3d4d_2p3p_3d4d/ints_h1.npy")
h2 = npzread("/Users/ayush/workspace/cMF_data/bimetallics/cr2_gagliardi/38_3d4d_2p3p_3d4d/ints_h2.npy")
ints = InCoreInts(h0, h1, h2)

clusters =  [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [11, 12, 13, 14, 15, 16, 17, 18, 19, 20], [21, 22, 23, 24, 25, 26], [27, 28, 29, 30, 31, 32], [33, 34, 35, 36, 37, 38]]
init_fspace =  [(3, 0), (3, 0), (3, 3), (3, 3), (3, 3)]

clusters = [MOCluster(i,collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)
println(n_orb(ints))
rdm1 = RDM1(n_orb(ints))
ansatze=[FCIAnsatz(length(ci),init_fspace[ci.idx][1],init_fspace[ci.idx][2]) for ci in clusters]
display(ansatze)
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_newton(ints, clusters, init_fspace,ansatze,rdm1, maxiter_oo = 400,
                           tol_oo=1e-6, 
                           tol_d1=1e-9, 
                           tol_ci=1e-11,
                           verbose=4, 
                           zero_intra_rots = true,
                           sequential=true)

ints = orbital_rotation(ints, U)

@save "data_cmf_cr2_gagliardi_38.jld2" clusters init_fspace ints d1 e_cmf U 