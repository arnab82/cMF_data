using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf
using ActiveSpaceSolvers
#C = npzread("mo_coeffs.npy")
#h0 = npzread("ints_h0.npy")
#h1 = npzread("ints_h1.npy")
#h2 = npzread("ints_h2.npy")
h0 = npzread("/home/arnab22/cMF_data/bimetallics/fe2_morokuma/631g/ints_h0.npy")
h1 = npzread("/home/arnab22/cMF_data/bimetallics/fe2_morokuma/631g/ints_h1.npy")
h2 = npzread("/home/arnab22/cMF_data/bimetallics/fe2_morokuma/631g/ints_h2.npy")
ints = InCoreInts(h0, h1, h2)

#Pa = npzread("Pa.npy")
#Pb = npzread("Pb.npy")
#@printf(" Input energy:    %12.8f\n", compute_energy(ints, RDM1(Pa, Pb)))


clusters =  [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56], [57, 58, 59, 60, 61], [62, 63, 64, 65, 66, 67], [68, 69, 70, 71, 72], [73, 74, 75, 76, 77, 78], [79, 80, 81, 82, 83, 84], [85, 86, 87, 88, 89, 90], [91, 92, 93, 94, 95, 96], [97, 98, 99, 100, 101, 102], [103, 104, 105, 106, 107, 108], [109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190]]
init_fspace =  [(56, 56), (5, 0), (3, 3), (5, 0), (3, 3), (3, 3), (3, 3), (3, 3), (3, 3), (3, 3), (0, 0)]
clusters = [MOCluster(i, collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)

rdm1 = RDM1(n_orb(ints))


# # Do CMF
#total number of electrons=164
ansatze = [FCIAnsatz(56,56,56),FCIAnsatz(5, 5,0), FCIAnsatz(6,3,3),FCIAnsatz(5,5,0),FCIAnsatz(6,3,3),FCIAnsatz(6,3,3),FCIAnsatz(6,3,3),FCIAnsatz(6,3,3),FCIAnsatz(6,3,3),FCIAnsatz(6,3,3),FCIAnsatz(82,0,0)]
# @time e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace,ansatze, RDM1(n_orb(ints)),
#     verbose=4, gconv=1e-8,tol_d1=1e-9,tol_ci=1e-11,sequential=true, max_iter_oo=100)
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_newton(ints, clusters, init_fspace,ansatze,rdm1, maxiter_oo = 400,
                           tol_oo=1e-6, 
                           tol_d1=1e-9, 
                           tol_ci=1e-11,
                           verbose=1, 
			   step_trust_region=1.5,
                           zero_intra_rots = true,
                           sequential=true,
			   trust_region=true)
ints = orbital_rotation(ints, U)

@save "data_fe2_morokuma_full_52_new.jld2" clusters init_fspace ints d1 e_cmf U
