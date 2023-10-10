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

h0 = npzread("/Users/ayush/workspace/cMF_data/new_pi_conj_sys/hbc.py.1693373.scr/ints_h0_hbc.npy")
h1 = npzread("/Users/ayush/workspace/cMF_data/new_pi_conj_sys/hbc.py.1693373.scr/ints_h1_hbc.npy")
h2 = npzread("/Users/ayush/workspace/cMF_data/new_pi_conj_sys/hbc.py.1693373.scr/ints_h2_hbc.npy")

ints = InCoreInts(h0, h1, h2)

#Hexabenzo[a,d,g,j,m,p]coronene


na = 24
nb = 24

clusters = [[16,32,46,43,28,19,12,6,4,10],[22,35,47,44,29,17],[20,31,39,42,34,18,9,3,5,11],[21,30,38,41,25,13],[15,27,40,37,26,14,8,2,1,7],[23,33,45,48,36,24]]
init_fspace = [(5,5),(3,3),(5,5),(3,3),(5,5),(3,3)]


rdm1 = RDM1(n_orb(ints))

# define clusters
clusters = [MOCluster(i,collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)
# @time e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace, RDM1(n_orb(ints)),
#     verbose=0, gconv=1e-6,tol_d1=1e-8,tol_ci=1e-9,sequential=true, max_iter_oo=100)
ansatze = [FCIAnsatz(10, 5, 5),FCIAnsatz(6, 3, 3), FCIAnsatz(10,5,5),FCIAnsatz(6, 3, 3),FCIAnsatz(10, 5, 5), FCIAnsatz(6,3,3)]

@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace,ansatze, rdm1,
                                                                   maxiter_oo   = 500, 
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
                                                                   orb_hessian=false)
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

@save "data_cmf_pi_hbc.jld2" clusters init_fspace ints d1 e_cmf U 