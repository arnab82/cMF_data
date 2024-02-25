using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf


@load "data_fe2_morokuma_full_52_new.jld2"
M = 50
init_fspace =  [(56, 56), (5, 0), (3, 3), (0,5), (3, 3), (3, 3), (3, 3), (3, 3), (3, 3), (3, 3), (0, 0)]
init_fspace = FermiCG.FockConfig(init_fspace)
cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [1,3,3,3,3,3,3,3,3,3,1], init_fspace, max_roots=M, verbose=1);

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=6


ci_vector = FermiCG.TPSCIstate(clusters, init_fspace, R=nroots)

ci_vector = FermiCG.add_spin_focksectors(ci_vector)
@save "cmfstate.jld2" ci_vector
eci, v = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);
@save "tpsstate.jld2" eci v
e2a = FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham,thresh_foi=1e-6)
e0a, v0a = FermiCG.tpsci_ci(v, cluster_ops, clustered_ham, incremental=true,
                            max_iter=1,
			    thresh_cipsi = 1e-2, 
                            thresh_spin  = 1e-3,
                            thresh_foi   = 1e-5,
                            max_mem_ci   = 200);
e2a_ = FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham,thresh_foi=1e-6)
@save "tpsstate.jld2" e0a v0a e2a_
