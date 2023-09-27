using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf
using BenchmarkTools
using ActiveSpaceSolvers
C = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/first_set/benzene_10mer/mo_coeffs.npy")
h0 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/first_set/benzene_10mer/ints_h0.npy")
h1 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/first_set/benzene_10mer/ints_h1.npy")
h2 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/first_set/benzene_10mer/ints_h2.npy")
ints = InCoreInts(h0, h1, h2)

Pa = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/first_set/benzene_10mer/Pa.npy")
Pb = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/first_set/benzene_10mer/Pb.npy")
@printf(" Input energy:    %12.8f\n", compute_energy(ints, RDM1(Pa, Pb)))

init_fspace =  [(3, 3) for i in 1:10]
clusters    =  [collect((i-1)*6+1:(i-1)*6+6) for i in 1:10]

clusters = [MOCluster(i, collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)

rdm1 = RDM1(n_orb(ints))


@time e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace, RDM1(n_orb(ints)),
    verbose=0, gconv=1e-8,tol_d1=1e-9,tol_ci=1e-11,sequential=true, max_iter_oo=100)
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace,RDM1(n_orb(ints)),
    verbose=0, gconv=1e-8,tol_d1=1e-9,tol_ci=1e-11,method="cg",sequential=true, max_iter_oo=100)

@time e_cmf, Ugd, d1gd = ClusterMeanField.cmf_oo_gd(ints, clusters, init_fspace, rdm1, 
                            maxiter_oo = 800,
                            tol_oo=1e-8, 
                            tol_d1=1e-9, 
                            tol_ci=1e-11, 
                            alpha=.1, 
                            sequential=true)

ansatze = [FCIAnsatz(6, 3, 3),FCIAnsatz(6, 3, 3), FCIAnsatz(6,3,3),FCIAnsatz(6, 3, 3),FCIAnsatz(6, 3, 3), FCIAnsatz(6,3,3),FCIAnsatz(6, 3, 3),FCIAnsatz(6, 3, 3), FCIAnsatz(6,3,3),FCIAnsatz(6,3,3)]
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

ints = orbital_rotation(ints, U)

@save "data_cmf_benzene_10mer.jld2" clusters init_fspace ints d1 e_cmf U 
