using RDM
using ClusterMeanField
using ActiveSpaceSolvers
using NPZ
using JLD2
using QCBase
using InCoreIntegrals

@load "/home/arnabbachhar/workspace/project_hessian/new_cmf/ClusterMeanField.jl/examples/tetracene_tetramer/ints_C_TT.jld2"

# define clusters
cluster_list = [collect(1:10), collect(11:20), collect(21:30), collect(31:40)]
clusters = [MOCluster(i,collect(cluster_list[i])) for i = 1:length(cluster_list)]
init_fspace = [ (5,5) for i in 1:4]
display(clusters)
d1=RDM1(n_orb(ints_sorted))
norb=n_orb(ints_sorted)
ansatze=[FCIAnsatz(10, 5, 5),FCIAnsatz(10, 5, 5), FCIAnsatz(10,5,5),FCIAnsatz(10, 5, 5)]

k=zeros(norb*(norb-1)÷2)
K = RDM.unpack_gradient(k, norb)
#K=reshape(k,(norb,norb)
Ui = exp(K)
        
tmp = RDM.orbital_rotation(ints_sorted,Ui)
ints_sorted.h1 .= tmp.h1
ints_sorted.h2 .= tmp.h2

tmp = RDM.orbital_rotation(d1,Ui)
d1.a .= tmp.a
d1.b .= tmp.b
println("*******************fci calculation**************************")
e, rdm1_dict, rdm2_dict = cmf_ci(ints_sorted, clusters, init_fspace,ansatze,d1, 
                                         maxiter_d1 = 500, 
                                         maxiter_ci = 500, 
                                         tol_d1 = 1e-7, 
                                         tol_ci = 1e-8, 
                                         verbose    = 4, 
                                         sequential = false)
println(e)
gd1, gd2 = ClusterMeanField.assemble_full_rdm(clusters, rdm1_dict, rdm2_dict)
println("*********************gradient******************************")
g_i = RDM.build_orbital_gradient(ints_sorted, gd1, gd2)
display(g_i)
println("**************analytical hessian*******************************")
analytical_hessian=RDM.build_orbital_hessian(ints_sorted,gd1,gd2)

proj_vec=ClusterMeanField.projection_vector(ansatze, clusters,norb)
anly_hessian=proj_vec'*analytical_hessian*proj_vec
display(anly_hessian)
println("*********************finite difference hessian*****************")
finite_hessian=ClusterMeanField.orbital_hessian_numerical(ints_sorted, clusters, k, init_fspace, d1)
fi_hessian=proj_vec'*finite_hessian*proj_vec
display(fi_hessian)
finite_hessian2=ClusterMeanField.orbital_hessian_finite_difference(ints_sorted, clusters, k, fspace, d1, ci_conv = 1e-8, verbose = 0,stepsize = 1e-5)
fi_hessian2=proj_vec'*finite_hessian2*proj_vec
display(fi_hessian2)
display(norm(fi_hessian-anly_hessian))
display(norm(fi_hessian2-anly_hessian))

@save "data_hessian_tetracene_tetramer.jld2"  anly_hessian fi_hessian fi_hessian2 proj_vec clusters init_fspace ints d1 