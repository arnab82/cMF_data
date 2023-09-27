using ClusterMeanField
using RDM
using InCoreIntegrals
using PyCall
using QCBase
using LinearAlgebra
using Printf
using ActiveSpaceSolvers

pyscf = pyimport("pyscf");
fcidump = pyimport("pyscf.tools.fcidump");
ctx = fcidump.read("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/PAHs/excited_tpsci_data/P2/biphenyl-2mer-1d-fcidump");
h = ctx["H1"];
g = ctx["H2"];
ecore = ctx["ECORE"];
g = pyscf.ao2mo.restore("1", g, size(h,2))
ints = InCoreInts(ecore,h,g)

rdm1 = RDM1(n_orb(ints))


na = 12
nb = 12

clusters_in    = [(1:6),(7:12),(13:18),(19:24)]
init_fspace = [(3,3),(3,3),(3,3),(3,3)]



# define clusters
clusters = [MOCluster(i,collect(clusters_in[i])) for i = 1:length(clusters_in)]
display(clusters)

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

ansatze = [FCIAnsatz(6, 3, 3),FCIAnsatz(6, 3, 3), FCIAnsatz(6,3,3),FCIAnsatz(6, 3, 3)]
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


@save "data_cmf_P2.jld2" clusters init_fspace ints d1 e_cmf U 












# max_roots = 100
# # Build Cluster basis
# cluster_bases = FermiCG.compute_cluster_eigenbasis(ints, clusters, verbose=0, max_roots=max_roots,
#         init_fspace=init_fspace, rdm1a=Da, rdm1b=Db);

# #cluster_bases = FermiCG.compute_cluster_est_basis(ints, clusters, Da, Db, thresh_schmidt=5e-5, init_fspace=init_fspace)

# #
# # Build ClusteredOperator
# clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

# # Build Cluster Operators
# cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

# # Add cmf hamiltonians for doing MP-style PT2 
# FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, Da, Db, verbose=0);



# nroots = 5

# ref_fock = FermiCG.FockConfig(init_fspace)
# ci_vector = FermiCG.ClusteredState(clusters, ref_fock, R=nroots)
# ci_vector[ref_fock][ClusterConfig([2,1,1,1])] = [0,1,0,0,0]
# ci_vector[ref_fock][ClusterConfig([1,2,1,1])] = [0,0,1,0,0]
# ci_vector[ref_fock][ClusterConfig([1,1,2,1])] = [0,0,0,1,0]
# ci_vector[ref_fock][ClusterConfig([1,1,1,2])] = [0,0,0,0,1]

# thresh_list = [0.005,0.002,0.001,0.0007,0.0005,0.0003,0.0001]


# for thresh_cipsi in thresh_list
#     e0, v0 = FermiCG.tpsci_ci(ci_vector, cluster_ops, clustered_ham,
#     			    thresh_cipsi=thresh_cipsi, # Threshold for adding to P-space
#     			    thresh_foi=1e-5,    # Threshold for keeping terms when defining FOIS    
#     			    thresh_asci=0.00,     # Threshold of P-space configs to search from
#     			    max_iter=10,
#     			    matvec=3);

#     e2a = FermiCG.compute_batched_pt2(v0, cluster_ops, clustered_ham, thresh_foi=1e-6)

#     println()
#     println("	*======TPSCI results======*")
#     @printf("TCI Thresh: %8.6f  Dim:%8d\n",thresh_cipsi,size(v0)[1])
#     println()
#     @printf("TCI %5s %12s %12s\n", "Root", "E(0)", "E(2)") 
#     for r in 1:nroots
#         @printf("TCI %5s %12.8f %12.8f\n",r, e0[r] + ecore, e0[r] + e2a[r] + ecore)
#     end
    
#     #global ci_vector = v0


#     for r in 1:nroots
#         @printf("TPSCI %5s %12.8f %12.8f\n",r, e0[r] + ecore, e0[r] + e2a[r] + ecore)
# 	display(v0,root=r)
#     end
# end

