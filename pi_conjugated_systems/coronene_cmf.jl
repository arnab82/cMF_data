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

h0 = npzread("/Users/ayush/workspace/fermicluster/ints_h0_coronene.npy")
h1 = npzread("/Users/ayush/workspace/fermicluster/ints_h1_coronene.npy")
h2 = npzread("/Users/ayush/workspace/fermicluster/ints_h2_coronene.npy")
# C= npzread("/Users/ayush/workspace/fermicluster/mo_coeffs_coronene.npy")

println(h0)
println(size(h1)[1])
println(size(h2)[1])



na = 12
nb = 12
# clusters = [[1,2,3,18,19,20],[4,5],[6,7,8,9,21,22],[10,11],[12,13,14,15,23,24],[16,17]]
clusters=[[1,2,3,4,5,6],[7,8],[9,10,11,12,13,14],[15,16],[17,18,19,20,21,22],[23,24]]
init_fspace = [(3,3),(1,1),(3,3),(1,1),(3,3),(1,1)]

ints = InCoreInts(h0, h1, h2)

rdm1 = RDM1(n_orb(ints))

# define clusters
clusters = [MOCluster(i,collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)

# @time e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace, RDM1(n_orb(ints)),
#     verbose=0, gconv=1e-8,tol_d1=1e-9,tol_ci=1e-11,sequential=true, max_iter_oo=100)
# @time e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace,RDM1(n_orb(ints)),
#     verbose=0, gconv=1e-8,tol_d1=1e-9,tol_ci=1e-11,method="cg",sequential=true, max_iter_oo=100)

# @time e_cmf, Ugd, d1gd = ClusterMeanField.cmf_oo_gd(ints, clusters, init_fspace, rdm1, 
#                             maxiter_oo = 800,
#                             tol_oo=1e-8, 
#                             tol_d1=1e-9, 
#                             tol_ci=1e-11, 
#                             alpha=.1, 
#                             sequential=true)
# @time e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace,rdm1,
#                                                                    maxiter_oo   = 25, 
#                                                                    maxiter_ci   = 200, 
#                                                                    maxiter_d1   = 200, 
#                                                                    verbose      = 0, 
#                                                                    tol_oo       = 1e-6, 
#                                                                    tol_d1       = 1e-8, 
#                                                                    tol_ci       = 1e-9, 
#                                                                    sequential   = true, 
#                                                                    alpha        = .1,
#                                                                    diis_start   = 1,
#                                                                    max_ss_size  = 24,
#                                                                    zero_intra_rots=true)

ansatze = [FCIAnsatz(6, 3, 3),FCIAnsatz(2,1,1), FCIAnsatz(6,3,3),FCIAnsatz(2,1,1),FCIAnsatz(6,3,3),FCIAnsatz(2,1,1)]
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace,ansatze, rdm1,
                                                                   maxiter_oo   = 250, 
                                                                   maxiter_ci   = 200, 
                                                                   maxiter_d1   = 200, 
                                                                   verbose      = 0, 
                                                                   tol_oo       = 1e-8, 
                                                                   tol_d1       = 1e-9, 
                                                                   tol_ci       = 1e-11, 
                                                                   sequential   = true, 
                                                                   diis_start   = 1,
                                                                   max_ss_size  = 24,
                                                                   orb_hessian=true,
                                                                   zero_intra_rots=false)
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace,ansatze, rdm1,
                                                                   maxiter_oo   = 250, 
                                                                   maxiter_ci   = 200, 
                                                                   maxiter_d1   = 200, 
                                                                   verbose      = 0, 
                                                                   tol_oo       = 1e-8, 
                                                                   tol_d1       = 1e-9, 
                                                                   tol_ci       = 1e-11, 
                                                                   sequential   = true, 
                                                                   alpha        = .1,
                                                                   diis_start   = 1,
                                                                   max_ss_size  = 24,
                                                                   orb_hessian=false,
                                                                   zero_intra_rots=true)
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace,ansatze, rdm1,
                                                                   maxiter_oo   = 250, 
                                                                   maxiter_ci   = 200, 
                                                                   maxiter_d1   = 200, 
                                                                   verbose      = 0, 
                                                                   tol_oo       = 1e-8, 
                                                                   tol_d1       = 1e-9, 
                                                                   tol_ci       = 1e-11, 
                                                                   sequential   = true, 
                                                                   diis_start   = 1,
                                                                   max_ss_size  = 24,
                                                                   orb_hessian=false,
                                                                   zero_intra_rots=false)
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace,ansatze, rdm1,
                                                                   maxiter_oo   = 500, 
                                                                   maxiter_ci   = 200, 
                                                                   maxiter_d1   = 200, 
                                                                   verbose      = 0, 
                                                                   tol_oo       = 1e-8, 
                                                                   tol_d1       = 1e-9, 
                                                                   tol_ci       = 1e-11, 
                                                                   sequential   = true, 
                                                                   diis_start   = 1,
                                                                   max_ss_size  = 24,
                                                                #    use_pyscf=false,
                                                                   orb_hessian=true,
                                                                   zero_intra_rots=true)


@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_newton(ints, clusters, init_fspace,ansatze,rdm1, maxiter_oo = 250,
                           tol_oo=1e-8, 
                           tol_d1=1e-9, 
                           tol_ci=1e-11,
                           verbose=4, 
                           zero_intra_rots = true,
                           sequential=true)
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_newton(ints, clusters, init_fspace,ansatze,rdm1, maxiter_oo = 250,
                           tol_oo=1e-8, 
                           tol_d1=1e-9, 
                           tol_ci=1e-11,
                           verbose=4, 
                           zero_intra_rots = false,
                           sequential=true)

ints = orbital_rotation(ints, U)

@save "data_cmf_coronene.jld2" clusters init_fspace ints d1 e_cmf U 




