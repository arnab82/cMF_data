using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf
using ActiveSpaceSolvers
C = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/ras_two_benzene/he_ras/mo_coeffs.npy")
h0 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/ras_two_benzene/he_ras/ints_h0.npy")
h1 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/ras_two_benzene/he_ras/ints_h1.npy")
h2 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/ras_two_benzene/he_ras/ints_h2.npy")
ints = InCoreInts(h0, h1, h2)

#Define clusters and intial Fock space for inital CMF calc for 2 orbs each He
clusters_in = [(1:8),(9:16)]
init_fspace = [(4,4),(4,4)]
rdm1 = zeros(size(ints.h1))
na=8
nb=8
#Define clusters now using ClusterMeanField code
clusters = [MOCluster(i,collect(clusters_in[i])) for i = 1:length(clusters_in)]
display(clusters)

print(size(ints.h1))

# ansatze = [RASCIAnsatz(8, 4,4, (2,4,2), max_h=4, max_p=4), RASCIAnsatz(8,4,4,(2,4,2), max_h=4, max_p=4)] #FCI type RASCI calculation
# ansatze = [RASCIAnsatz(8,4,4, (2,4,2), max_h=0, max_p=0), RASCIAnsatz(8,4,4,(2,4,2), max_h=0, max_p=0)] #CASCI type RASCI calculation
# ansatze = [RASCIAnsatz(8, 4,4, (2,4,2), max_h=2, max_p=2), RASCIAnsatz(8,4,4,(2,4,2), max_h=2, max_p=2)] # double excitationRASCI calculation
ansatze = [RASCIAnsatz(8, 4,4, (2,4,2), max_h=3, max_p=3), RASCIAnsatz(8,4,4,(2,4,2), max_h=3, max_p=3)] # excitation RASCI calculation
# ansatze = [RASCIAnsatz(8, 4,4, (2,4,2), max_h=1, max_p=1), RASCIAnsatz(8,4,4,(2,4,2), max_h=1, max_p=1)] #single excitations type RASCI calculation
display(ansatze)

rdm1 = RDM1(n_orb(ints))
e_cmf, U_cmf, d1  = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, ansatze,  rdm1, maxiter_oo = 400,
                                                 tol_oo=1e-7,tol_d1=1e-9, tol_ci=1e-11, verbose=0, diis_start=1,sequential=true);
ints_cmf = orbital_rotation(ints,U_cmf)

e_cmf, U_n, d1_n = ClusterMeanField.cmf_oo_newton2(ints, clusters, init_fspace, ansatze,rdm1, maxiter_oo = 2500,
                           tol_oo=1e-7, 
                           tol_d1=1e-9, 
                           tol_ci=1e-11,
                           verbose=4, 
                           zero_intra_rots = true,
                           sequential=true)
ints_cmf = orbital_rotation(ints,U_n)

