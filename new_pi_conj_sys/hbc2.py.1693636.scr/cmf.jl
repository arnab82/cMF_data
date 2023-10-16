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

h0 = npzread("/Users/ayush/workspace/cMF_data/new_pi_conj_sys/hbc2.py.1693636.scr/ints_h0_hbc2.npy")
h1 = npzread("/Users/ayush/workspace/cMF_data/new_pi_conj_sys/hbc2.py.1693636.scr/ints_h1_hbc2.npy")
h2 = npzread("/Users/ayush/workspace/cMF_data/new_pi_conj_sys/hbc2.py.1693636.scr/ints_h2_hbc2.npy")

ints = InCoreInts(h0, h1, h2)

#Hexabenzo[a,d,g,j,m,p]coronene


na = 21
nb = 21

clusters = [[7,15,27,37,25,13],[19,31,39,29,17,9],[23,35,41,33,21,11],[22,34,42,36,24,12],[18,30,40,32,20,10],[14,26,38,28,16,8],[2,4,6,5,3,1]]
init_fspace = [(3,3),(3,3),(3,3),(3,3),(3,3),(3,3),(3,3)]


rdm1 = RDM1(n_orb(ints))

# define clusters
clusters = [MOCluster(i,collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)
# @time e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace, RDM1(n_orb(ints)),
#     verbose=0, gconv=1e-6,tol_d1=1e-8,tol_ci=1e-9,sequential=true, max_iter_oo=100)
ansatze = [FCIAnsatz(6, 3, 3),FCIAnsatz(6, 3, 3), FCIAnsatz(6,3,3),FCIAnsatz(6, 3, 3),FCIAnsatz(6, 3, 3), FCIAnsatz(6,3,3),FCIAnsatz(6,3,3)]

# @time e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace,ansatze, rdm1,
#                                                                    maxiter_oo   = 500, 
#                                                                    maxiter_ci   = 200, 
#                                                                    maxiter_d1   = 200, 
#                                                                    verbose      = 0, 
#                                                                    tol_oo       = 1e-8, 
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
                                                                   tol_oo       = 1e-8, 
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
                                                                   tol_oo       = 1e-8, 
                                                                   tol_d1       = 1e-9, 
                                                                   tol_ci       = 1e-11, 
                                                                   sequential   = true, 
                                                                   diis_start   = 1,
                                                                   max_ss_size  = 24,
                                                                #    use_pyscf=false,
                                                                   orb_hessian=true,
                                                                   zero_intra_rots=true)


@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_newton(ints, clusters, init_fspace,ansatze,rdm1, maxiter_oo = 400,
                           tol_oo=1e-8, 
                           tol_d1=1e-9, 
                           tol_ci=1e-11,
                           verbose=4, 
                           zero_intra_rots = true,
                           sequential=true)
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_newton(ints, clusters, init_fspace,ansatze,rdm1, maxiter_oo = 400,
                           tol_oo=1e-8, 
                           tol_d1=1e-9, 
                           tol_ci=1e-11,
                           verbose=4, 
                           zero_intra_rots = false,
                           sequential=true)


ints = orbital_rotation(ints, U)

@save "data_cmf_pi_hbc2.jld2" clusters init_fspace ints d1 e_cmf U 