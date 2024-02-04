using QCBase
using ClusterMeanField 
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using ActiveSpaceSolvers

h0 = npzread("/Users/ayush/workspace/cMF_data/bimetallics/cr2_gagliardi/3d4d_2p3p_3d4d_docc_virt/ints_h0.npy")
h1 = npzread("/Users/ayush/workspace/cMF_data/bimetallics/cr2_gagliardi/3d4d_2p3p_3d4d_docc_virt/ints_h1.npy")
h2 = npzread("/Users/ayush/workspace/cMF_data/bimetallics/cr2_gagliardi/3d4d_2p3p_3d4d_docc_virt/ints_h2.npy")
ints = InCoreInts(h0, h1, h2)
init_fspace =  [(3, 0), (3, 0), (4, 4), (4, 4), (4, 4), (14, 14), (0, 0)]
clusters    =  [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [11, 12, 13, 14, 15, 16, 17, 18, 19, 20], [21, 22, 23, 24, 25, 26, 27], [28, 29, 30, 31, 32, 33, 34], [35, 36, 37, 38, 39, 40, 41], [42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55], [56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75]]

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
                           verbose=0,
                           step_trust_region=4.5 ,
                           zero_intra_rots = true,
                           sequential=true,
                           trust_region=true)

ints = orbital_rotation(ints, U)

@save "data_cmf_cr2_gagliardi_75_3d4d.jld2" clusters init_fspace ints d1 e_cmf U 