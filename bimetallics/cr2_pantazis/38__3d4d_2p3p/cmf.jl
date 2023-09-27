using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf

C = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/bimetallics/cr2_pantazis/38__3d4d_2p3p/mo_coeffs.npy")
h0 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/bimetallics/cr2_pantazis/38__3d4d_2p3p/ints_h0.npy")
h1 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/bimetallics/cr2_pantazis/38__3d4d_2p3p/ints_h1.npy")
h2 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/bimetallics/cr2_pantazis/38__3d4d_2p3p/ints_h2.npy")
ints = InCoreInts(h0, h1, h2)

Pa = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/bimetallics/cr2_pantazis/38__3d4d_2p3p/Pa.npy")
Pb = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/bimetallics/cr2_pantazis/38__3d4d_2p3p/Pb.npy")
@printf(" Input energy:    %12.8f\n", compute_energy(ints, RDM1(Pa, Pb)))

 init_fspace =  [(5, 2), (2, 5), (3, 3), (3, 3), (3, 3)]
 clusters    =  [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [11, 12, 13, 14, 15, 16, 17, 18, 19, 20], [21, 22, 23, 24, 25, 26], [27, 28, 29, 30, 31, 32], [33, 34, 35, 36, 37, 38]]

clusters = [MOCluster(i, collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)

d1 = RDM1(n_orb(ints))


rdm1 = RDM1(n_orb(ints))

ansatze = [FCIAnsatz(10,5,2), FCIAnsatz(10,2,5),FCIAnsatz(6,3,3),FCIAnsatz(6,3,3),FCIAnsatz(6,3,3)]
# @time e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace, RDM1(n_orb(ints)),
#     verbose=0, gconv=1e-8,tol_d1=1e-9,tol_ci=1e-11,sequential=true, max_iter_oo=100)
# @time e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace,RDM1(n_orb(ints)),
#     verbose=0, gconv=1e-8,tol_d1=1e-9,tol_ci=1e-11,method="cg",sequential=true, max_iter_oo=100)

# @time e_cmf, Ugd, d1gd = ClusterMeanField.cmf_oo_gd(ints, clusters, init_fspace, rdm1, 
#                             maxiter_oo = 1800,
#                             tol_oo=1e-8, 
#                             tol_d1=1e-9, 
#                             tol_ci=1e-11, 
#                             alpha=.1, 
#                             sequential=true)
# @time e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace,ansatze, rdm1,
#                                                                    maxiter_oo   = 500, 
#                                                                    maxiter_ci   = 200, 
#                                                                    maxiter_d1   = 200, 
#                                                                    verbose      = 0, 
#                                                                    tol_oo       = 1e-7, 
#                                                                    tol_d1       = 1e-9, 
#                                                                    tol_ci       = 1e-11, 
#                                                                    sequential   = true, 
#                                                                    alpha        = .1,
#                                                                    diis_start   = 1,
#                                                                    max_ss_size  = 24,
#                                                                    orb_hessian=false)
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace,ansatze, rdm1,
                                                                   maxiter_oo   = 500, 
                                                                   maxiter_ci   = 200, 
                                                                   maxiter_d1   = 200, 
                                                                   verbose      = 0, 
                                                                   tol_oo       = 1e-7, 
                                                                   tol_d1       = 1e-9, 
                                                                   tol_ci       = 1e-11, 
                                                                   sequential   = true, 
                                                                   diis_start   = 1,
                                                                   max_ss_size  = 24,
                                                                   orb_hessian=true,
                                                                   zero_intra_rots=false)
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace,ansatze, rdm1,
                                                                   maxiter_oo   = 500, 
                                                                   maxiter_ci   = 200, 
                                                                   maxiter_d1   = 200, 
                                                                   verbose      = 0, 
                                                                   tol_oo       = 1e-7, 
                                                                   tol_d1       = 1e-9, 
                                                                   tol_ci       = 1e-11, 
                                                                   sequential   = true, 
                                                                   diis_start   = 1,
                                                                   max_ss_size  = 24,
                                                                #    use_pyscf=false,
                                                                   orb_hessian=true,
                                                                   zero_intra_rots=true)


@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_newton2(ints, clusters, init_fspace,ansatze,rdm1, maxiter_oo = 400,
                           tol_oo=1e-7, 
                           tol_d1=1e-9, 
                           tol_ci=1e-11,
                           verbose=4, 
                           zero_intra_rots = true,
                           sequential=true)
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_newton2(ints, clusters, init_fspace,ansatze,rdm1, maxiter_oo = 400,
                           tol_oo=1e-7, 
                           tol_d1=1e-9, 
                           tol_ci=1e-11,
                           verbose=4, 
                           zero_intra_rots = false,
                           sequential=true)
ints = orbital_rotation(ints, U)
@save "cr2_pantazis_38__3d4d_2p3p_data_cmf.jld2" clusters init_fspace ints d1 e_cmf U 
