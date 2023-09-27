using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf
using BenchmarkTools
h0 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/first_set/tetraphenyl_ethylene/tetraphenyl-ethylene_scf_integrals_h0.npz.npy")
h1 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/first_set/tetraphenyl_ethylene/tetraphenyl-ethylene_scf_integrals_h1.npz.npy")
h2 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/first_set/tetraphenyl_ethylene/tetraphenyl-ethylene_scf_integrals_h2.npz.npy");
ints = InCoreInts(h0, h1, h2)
print(size(h0), size(h1), size(h2))
# Define clusters

clusters_in = [
    (1:6),   # Benzene 1
    (7:12),  # Benzene 2
    (13:18), # Benzene 3
    (19:24), # Benzene 4
    (25:26)  # ethylene
]
clusters = [MOCluster(i,collect(clusters_in[i])) for i = 1:length(clusters_in)]

init_fspace = [
    (3,3),
    (3,3),
    (3,3),
    (3,3),
    (1,1)
];
print(clusters)
print(init_fspace)

rdm1 = RDM1(n_orb(ints))
# # # Do CMF
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

ansatze = [FCIAnsatz(6, 3, 3),FCIAnsatz(6, 3, 3), FCIAnsatz(6,3,3),FCIAnsatz(6, 3, 3),FCIAnsatz(2,1,1)]
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


@save "data_cmf_tetraphenyl_ethylene.jld2" clusters init_fspace ints d1 e_cmf U 
