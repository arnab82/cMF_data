using QCBase
using ClusterMeanField
using ActiveSpaceSolvers
using RDM
using InCoreIntegrals
using PyCall
using JLD2
using NPZ
using Printf

molecule = "
C       -6.8604981940     -0.4376683883      0.0000000000
C       -7.0417938838      0.9473143671      0.0000000000
C       -5.5765818320     -0.9841959654      0.0000000000
C       -5.9160648559      1.7742824642      0.0000000000
C       -4.4255614400     -0.1700917199      0.0000000000
C       -4.6336867757      1.2242295397      0.0000000000
C       -3.0514687259     -0.7633813258      0.0000000000
C       -2.8520613255     -2.1570405955      0.0000000000
C       -1.5640941394     -2.6887327507      0.0000000000
C       -0.4451033014     -1.8584342989      0.0000000000
C       -0.5927492072     -0.4581981628      0.0000000000
C       -1.9041039939      0.0505338987      0.0000000000
H       -8.0463585128      1.3756739224      0.0000000000
H       -5.4821758879     -2.0696091478      0.0000000000
H       -6.0325833442      2.8604965477      0.0000000000
H       -3.7849879120      1.9076345469      0.0000000000
H       -3.6983932318     -2.8425398416      0.0000000000
H       -1.4298106521     -3.7726953568      0.0000000000
H        0.5426199813     -2.3168156796      0.0000000000
H       -2.0368540411      1.1280906363      0.0000000000
H       -7.7260669052     -1.1042415323      0.0000000000
H        0.22700           0.16873           0.00000
"
atoms = []
for (li,line) in enumerate(split(rstrip(lstrip(molecule)), "\n"))
    l = split(line)
    push!(atoms, Atom(li, l[1], parse.(Float64,l[2:4])))
end

basis = "cc-pvdz"
# Create FermiCG.Molecule type
pyscf = pyimport("pyscf")
tools = pyimport("pyscf.tools")
fcidump = pyimport("pyscf.tools.fcidump");

pymol = pyscf.gto.Mole(atom=molecule,
                       symmetry = false, spin =0,charge=0,
                       basis = basis)
pymol.build()

h0 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/ras_two_benzene/ints_h0.npy")
h1 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/ras_two_benzene/ints_h1.npy")
h2 = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/ras_two_benzene/ints_h2.npy")
ints = InCoreInts(h0, h1, h2)

Cact = npzread("/Users/ayush/workspace/cmf/project_hessian_debug/ClusterMeanField.jl/examples/ras_two_benzene/mo_coeffs.npy")

clusters_in    = [(1:6),(7:12)]
n_clusters = 2
cluster_list = [collect(1:6), collect(7:12)]
clusters = [MOCluster(i,collect(cluster_list[i])) for i = 1:length(cluster_list)]
init_fspace = [ (3,3) for i in 1:n_clusters]
display(clusters)
# ansatze = [RASCIAnsatz(6, 3, 3, (1,4,1), max_h=2, max_p=2), RASCIAnsatz(6,3,3,(1,4,1), max_h=2, max_p=2)] #FCI type RASCI calculation
ansatze = [RASCIAnsatz(6, 3, 3, (1,4,1), max_h=1, max_p=1), RASCIAnsatz(6,3,3,(1,4,1), max_h=1, max_p=1)] #Single excitation RASCI calculation
# ansatze = [RASCIAnsatz(6, 3, 3, (1,4,1), max_h=0, max_p=0), RASCIAnsatz(6,3,3,(1,4,1), max_h=0, max_p=0)] #CASCI type RASCI calculation
display(ansatze)

rdm1 = zeros(size(ints.h1))
# run cmf_oo
# e_cmf, U_cmf, d1  = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, ansatze, RDM1(rdm1, rdm1), maxiter_oo = 400, tol_oo=1e-8,tol_d1=1e-9, tol_ci=1e-10, verbose=0, diis_start=1);
# ints_cmf = orbital_rotation(ints,U_cmf)
# C_cmf = Cact*U_cmf
e_cmf, U_n, d1_n = ClusterMeanField.cmf_oo_newton2(ints, clusters, init_fspace, ansatze,RDM1(rdm1, rdm1), maxiter_oo = 2500,
                           tol_oo=1e-8, 
                           tol_d1=1e-9, 
                           tol_ci=1e-10,
                           verbose=4, 
                           zero_intra_rots = true,
                           sequential=true)
ints_cmf = orbital_rotation(ints,U_n)
C_cmf = Cact*U_n













"""

@testset "CMF_RAS" begin
    atoms = []
    push!(atoms,Atom(1,"He",[-2,0,0]))
    push!(atoms,Atom(2,"He",[2,0,0]))
    push!(atoms,Atom(3,"He",[4,0,0]))
    push!(atoms,Atom(4,"He",[6,0,0]))
    basis="cc-pvdz"
    mol     = Molecule(0,1,atoms,basis)
    mf = ClusterMeanField.pyscf_do_scf(mol)
    nbas = size(mf.mo_coeff)[1]
    ints = ClusterMeanField.pyscf_build_ints(mol,mf.mo_coeff, zeros(nbas,nbas));
    @printf(" HF  Energy: %12.8f\n", mf.e_tot)
    C = mf.mo_coeff


    rdm_mf = mf.mo_coeff'*mf.get_ovlp()*mf.make_rdm1()*mf.get_ovlp()*mf.mo_coeff / 2
    rdm1 = RDM1(rdm_mf, rdm_mf)
    @printf(" Should be E(HF):  %12.8f\n", compute_energy(ints, rdm1, RDM2(rdm1)))
    Cl = ClusterMeanField.localize(mf.mo_coeff,"lowdin",mf)
    ClusterMeanField.pyscf_write_molden(mol,Cl,filename="lowdin_He4.molden")
    S = ClusterMeanField.get_ovlp(mf)
    U =  C' * S * Cl
    println(" Build Integrals")
    flush(stdout)
    ints = orbital_rotation(ints,U)
    rdm1 = orbital_rotation(rdm1,U)
    println(" done.")
    flush(stdout)

    init_fspace=[(2,2),(1,1)]
    cluster=[(1:10),(11:20)]
    clusters = [MOCluster(i,collect(cluster[i])) for i = 1:length(cluster)]
    display(clusters)
    ansatze=[RASCIAnsatz(10,2,2,(1,5,4)),RASCIAnsatz(10,2,2,(1,5,4))]

    e_cmf, U = cmf_oo(ints, clusters, init_fspace,ansatze, rdm1, 
                              verbose=0, gconv=1e-6, method="cg",sequential=true)
    @test isapprox(e_cmf, -8.660402773572, atol=1e-10)                       
    e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace,ansatze, rdm1,
                        maxiter_oo   = 500,
                        maxiter_ci   = 200,
                        maxiter_d1   = 200,
                        verbose      = 0,
                        tol_oo       = 1e-6,
                        tol_d1       = 1e-9,
                        tol_ci       = 1e-11,
                        sequential   = true,
                        diis_start   = 1,
                        max_ss_size  = 24)
    @test isapprox(e_cmf, -8.660402773572, atol=1e-10)   
end"""