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

h0 = npzread("/home/arnabbachhar/workspace/project_hessian/new_cmf/cMF_data/new_pi_conj_sys/hbc.py.1693373.scr/ints_h0_hbc.npy")
h1 = npzread("/home/arnabbachhar/workspace/project_hessian/new_cmf/cMF_data/new_pi_conj_sys/hbc.py.1693373.scr/ints_h1_hbc.npy")
h2 = npzread("/home/arnabbachhar/workspace/project_hessian/new_cmf/cMF_data/new_pi_conj_sys/hbc.py.1693373.scr/ints_h2_hbc.npy")

ints = InCoreInts(h0, h1, h2)



na = 24
nb = 24

clusters = [[43,46,32,16,19,28],[35,22,17,29,44,47],[20,18,34,42,39,31],[30,38,41,25,13,21],[14,15,27,40,37,26],[24,23,33,45,48,36],[12,6,5,11,9,3,1,7,8,2,4,10]]
init_fspace = [(3,3),(3,3),(3,3),(3,3),(3,3),(3,3),(6,6)]


rdm1 = RDM1(n_orb(ints))

# define clusters
clusters = [MOCluster(i,collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace, RDM1(n_orb(ints)),
    verbose=0, gconv=1e-6,tol_d1=1e-8,tol_ci=1e-9,sequential=true, max_iter_oo=100)
ansatze = [FCIAnsatz(6, 3, 3),FCIAnsatz(6, 3, 3), FCIAnsatz(6,3,3),FCIAnsatz(6, 3, 3),FCIAnsatz(6, 3, 3), FCIAnsatz(6,3,3),FCIAnsatz(12,6,6)]

# @time e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace,ansatze, rdm1,
#                                                                    maxiter_oo   = 500, 
#                                                                    maxiter_ci   = 200, 
#                                                                    maxiter_d1   = 200, 
#                                                                    verbose      = 0, 
#                                                                    tol_oo       = 1e-6, 
#                                                                    tol_d1       = 1e-8, 
#                                                                    tol_ci       = 1e-9, 
#                                                                    sequential   = true, 
#                                                                    diis_start   = 1,
#                                                                    max_ss_size  = 24,
#                                                                    orb_hessian=true,
#                                                                    zero_intra_rots=false)

ints = orbital_rotation(ints, U)

@save "data_cmf_pi_hbc_7cluster.jld2" clusters init_fspace ints d1 e_cmf U 