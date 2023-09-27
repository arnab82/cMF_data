using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf
using ActiveSpaceSolvers

h0 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/noble_gas_clusters/helium/ints_h0_he08_08.npy")
h1 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/noble_gas_clusters/helium/ints_h1_he08_08.npy")
h2 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/noble_gas_clusters/helium/ints_h2_he08_08.npy")
ints = InCoreInts(h0, h1, h2)

#Define clusters and intial Fock space for inital CMF calc for 5 orbs each He
clusters_in = [(1:5),(6:10), (11:15), (16:20), (21:25), (26:30), (31:35), (36:40)]
init_fspace = [(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
rdm1 = zeros(size(ints.h1))
na=8
nb=8
#Define clusters now using FermiCG code
clusters = [MOCluster(i,collect(clusters_in[i])) for i = 1:length(clusters_in)]
display(clusters)

print(size(ints.h1))
rdm1 =RDM1(n_orb(ints))

#do a CMF calculation to optimize cluster orbitals

@time e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace, RDM1(n_orb(ints)),
    verbose=0, gconv=1e-8,tol_d1=1e-9,tol_ci=1e-11,sequential=true, max_iter_oo=100)
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace,RDM1(n_orb(ints)),
    verbose=0, gconv=1e-8,tol_d1=1e-9,tol_ci=1e-11,method="cg",sequential=true, max_iter_oo=100)

@time e_cmf, Ugd, d1gd = ClusterMeanField.cmf_oo_gd(ints, clusters, init_fspace, rdm1, 
                            maxiter_oo = 1900,
                            tol_oo=1e-8, 
                            tol_d1=1e-9, 
                            tol_ci=1e-11, 
                            alpha=.1, 
                            sequential=true)

ansatze = [ FCIAnsatz(5,1,1),FCIAnsatz(5,1,1),FCIAnsatz(5,1,1),FCIAnsatz(5,1,1),FCIAnsatz(5,1,1),FCIAnsatz(5,1,1), FCIAnsatz(5,1,1),FCIAnsatz(5,1,1)]
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


@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_newton2(ints, clusters, init_fspace,ansatze,rdm1, maxiter_oo = 400,
                           tol_oo=1e-8, 
                           tol_d1=1e-9, 
                           tol_ci=1e-11,
                           verbose=4, 
                           zero_intra_rots = true,
                           sequential=true)
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_newton2(ints, clusters, init_fspace,ansatze,rdm1, maxiter_oo = 400,
                           tol_oo=1e-8, 
                           tol_d1=1e-9, 
                           tol_ci=1e-11,
                           verbose=4, 
                           zero_intra_rots = false,
                           sequential=true)


cmf_ints=orbital_rotation(ints,U)

@save "data_cmf_he08_08.jld2" clusters init_fspace cmf_ints d1 e_cmf U 