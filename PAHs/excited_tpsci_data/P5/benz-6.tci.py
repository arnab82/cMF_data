import sys, os
import numpy as np
import scipy
import itertools
import time
from math import factorial
import copy as cp
import sys
import tools 

#will need to impmort fermicluster package
from fermicluster import *
from pyscf_helper import *
import pyscf
from pyscf import gto, scf, ao2mo, molden, lo, mo_mapping, mcscf
#if errors occur then import molden from tools and check version of pyscf

pyscf.lib.num_threads(1) #with degenerate states and multiple processors there can be issues
np.set_printoptions(suppress=True, precision=3, linewidth=1500)


molecule = '''
  C      4.96969861      2.86944105      0.00000000
  C      3.76923841      3.57603699      0.00000000
  C      2.53681227      2.89535727      0.00000000
  C      2.57923172      1.48927197      0.00000000
  C      3.77589675      0.74907202      0.00000000
  C      4.98146389      1.47660211      0.00000000
  C      3.76923841     -3.57603699      0.00000000
  C      4.96969861     -2.86944105      0.00000000
  C      2.53681227     -2.89535727      0.00000000
  C      4.98146389     -1.47660211      0.00000000
  C      2.57923172     -1.48927197      0.00000000
  C      3.77589675     -0.74907202      0.00000000
  C      0.00000000     -5.73784656      0.00000000
  C      1.21207936     -5.05143158      0.00000000
  C     -1.21207936     -5.05143158      0.00000000
  C      1.23908502     -3.64398497      0.00000000
  C     -1.23908502     -3.64398497      0.00000000
  C      0.00000000     -2.97739002      0.00000000
  C     -3.76923841     -3.57603699      0.00000000
  C     -2.53681227     -2.89535727      0.00000000
  C     -4.96969861     -2.86944105      0.00000000
  C     -2.57923172     -1.48927197      0.00000000
  C     -4.98146389     -1.47660211      0.00000000
  C     -3.77589675     -0.74907202      0.00000000
  C     -3.77589675      0.74907202      0.00000000
  C     -2.57923172      1.48927197      0.00000000
  C     -4.98146389      1.47660211      0.00000000
  C     -2.53681227      2.89535727      0.00000000
  C     -4.96969861      2.86944105      0.00000000
  C     -3.76923841      3.57603699      0.00000000
  C      0.00000000      2.97739002      0.00000000
  C      1.23908502      3.64398497      0.00000000
  C     -1.23908502      3.64398497      0.00000000
  C      1.21207936      5.05143158      0.00000000
  C     -1.21207936      5.05143158      0.00000000
  C      0.00000000      5.73784656      0.00000000
  H     -5.91553326      3.41562910      0.00000000
  H     -3.80787403      4.66429963      0.00000000
  H      0.00000000      1.89135768      0.00000000
  H      2.13496388      5.62940829      0.00000000
  H     -2.13496388      5.62940829      0.00000000
  H      0.00000000      6.83007222      0.00000000
  H      5.91553326      3.41562910      0.00000000
  H      3.80787403      4.66429963      0.00000000
  H      1.63845757      0.94703610      0.00000000
  H      5.94345659      0.96667697      0.00000000
  H     -1.63845757      0.94703610      0.00000000
  H     -5.94345659      0.96667697      0.00000000
  H      5.94345659     -0.96667697      0.00000000
  H      5.91553326     -3.41562910      0.00000000
  H      3.80787403     -4.66429963      0.00000000
  H      1.63845757     -0.94703610      0.00000000
  H      0.00000000     -6.83007222      0.00000000
  H      2.13496388     -5.62940829      0.00000000
  H     -2.13496388     -5.62940829      0.00000000
  H      0.00000000     -1.89135768      0.00000000
  H     -3.80787403     -4.66429963      0.00000000
  H     -5.91553326     -3.41562910      0.00000000
  H     -1.63845757     -0.94703610      0.00000000
  H     -5.94345659     -0.96667697      0.00000000
'''
cas_nel = 36
cas_norb = 36



#PYSCF inputs
mol = gto.Mole(atom=molecule,
    symmetry = True,basis = 'ccpvdz' )
mol.build()
print("symmertry: ",mol.topgroup)

#SCF
mf = scf.RHF(mol)
mf.verbose = 4
mf.conv_tol = 1e-12
mf.conv_tol_grad = 1e-9
mf.run(max_cycle=200)

## Active space selection


h,ecore,g,C = get_pi_space(mol,mf,cas_norb,cas_nel,local=True)
print("core %16.8f"%ecore)



mc = mulliken_ordering(mol,h.shape[0],C)
idx = np.where(mc>.9)[1]  #gives index map from atom to local orbital corresponding to that orbital

# Reorder
h,g = reorder_integrals(idx,h,g)
print(h)
C = C[:,idx] # make sure u reorder this too
molden.from_mo(mol, 'cas.molden', C)


from pyscf import tools
tools.fcidump.from_integrals('benz-6-fcidump', h, g, cas_norb,
                                 cas_nel, nuc=ecore, ms=0)


## TPSCI
blocks = [range(0,6),range(6,12),range(12,18),range(18,24),range(24,30),range(30,36)]
init_fspace = ((3,3),(3,3),(3,3),(3,3),(3,3),(3,3))



if 0:
    # Initialize the CMF solver. 
    oocmf = CmfSolver(h, g, ecore, blocks, init_fspace,C,max_roots=100,cs_solver=0) #cs_solver,0 for our FCI and 1 for pyscf FCI solver.
    oocmf.init() # runs a single step CMF calculation
    oocmf.optimize_orbitals()  # optimize the orbitals using gradient
    oocmf.form_extra_fspace()  #form excited fock space configurations

    clustered_ham = oocmf.clustered_ham  # clustered_ham used for TPSCI calculation
    ci_vector = oocmf.ci_vector   # lowest energy TPS using the reference Fock space
    h = oocmf.h
    g = oocmf.g
    C = oocmf.C



    ##
    # Excited State TPSCI for the lowest 5 roots. The guess is generated using a CIS type guess, hence there can be other CT states lower which might need different
    #    initialization
    ##
    nroots = 7

    # STEP 1: Expand the CIS space and generate roots
    ci_vector_s = ci_vector.copy()
    ci_vector_s.add_single_excitonic_states(clustered_ham.clusters)
    H = build_full_hamiltonian_parallel2(clustered_ham, ci_vector_s)
    e,v = np.linalg.eigh(H)
    idx = e.argsort()
    e = e[idx]
    v = v[:,idx]

    #STEP 2: Store the first n roots into a clustered_state and prune
    all_vecs = []
    for rn in range(nroots):
        vec = ci_vector_s.copy()
        vec.zero()
        vec.set_vector(v[:,rn])
        vec.clip(1e-5)
        print("Root:%4d     Energy:%12.8f   Gap:%12.8f  CI Dim: %4i "%(rn,e[rn].real,e[rn].real-e[0].real,len(vec)))
        #vec.print_configs()
        all_vecs.append(vec)

    #STEP 3: Combine all the vecs for each roots into 1 single ci_vector space.
    # Note the coefficients are meanning less here since its multiple states
    for vi,vec in enumerate(all_vecs):
        ci_vector.add(vec)
    ci_vector.zero()
    ci_vector.print()

    #Exta: Print excited state energies after pruning. Just to make sure we have not lost any state
    H = build_full_hamiltonian_parallel2(clustered_ham, ci_vector, nproc=None)
    e,v = np.linalg.eigh(H)
    idx = e.argsort()
    e = e[idx]
    v = v[:,idx]
    for rn in range(nroots):
        print("Root:%4d     Energy:%12.8f   Gap:%12.8f"%(rn,e[rn].real,e[rn].real-e[0].real))


    #STEP 4: Run the Excited State-TPSCI and analyze results
    time1 = time.time()
    ci_vector, pt_vector, e0, e2  = ex_tp_cipsi(ci_vector, clustered_ham,  
        thresh_cipsi    = 1e-5, 
        thresh_conv     = 1e-8, 
        max_iter        = 30, 
        n_roots         = nroots,
        thresh_asci     = 1e-2,
        nbody_limit     = 4, 
        pt_type         = 'en',
        thresh_search   = 2e-6, 
        shared_mem      = 50e9,
        batch_size      = 1,
        matvec          = 4,
        nproc           = None)
    time2 = time.time()

    for rn in range(nroots):
        print("Root:%4d     Var Energy:%12.8f   Gap:%12.8f  CI Dim: %4i "%(rn,e0[rn].real,e0[rn].real-e0[0].real,len(ci_vector)))
    for rn in range(nroots):
        print("Root:%4d     PT  Energy:%12.8f   Gap:%12.8f  CI Dim: %4i "%(rn,e2[rn],e2[rn]-e2[0],len(ci_vector)))

    print("Time spent in the Ex-TPSCI code%16.8f"%(time2-time1))


    for rn in range(nroots):

        e2,_ = compute_pt2_correction(ci_vector[rn], clustered_ham, e0[rn],
                thresh_asci     = 0,
                thresh_search   = 1e-6,
                batch_size      = 10,
                shared_mem          = 2e11,
                pt_type         = 'mp',
                nproc = 24,
                matvec          = 4)

        print("TCI   Root:%4d     Var Energy:%12.8f   Gap:%12.8f  CI Dim: %4i "%(rn,e0[rn].real+ecore,e0[rn].real-e0[0].real,len(ci_vector)))
        print("TCIPT Root:%4d     Var Energy:%12.8f "%(rn,e2+ecore))


    #STEP 4: Run the Excited State-TPSCI and analyze results
    time1 = time.time()
    ci_vector, pt_vector, e0, e2  = ex_tp_cipsi(ci_vector, clustered_ham,  
        thresh_cipsi    = 5e-6, 
        thresh_conv     = 1e-8, 
        max_iter        = 30, 
        n_roots         = nroots,
        thresh_asci     = 1e-2,
        nbody_limit     = 4, 
        pt_type         = 'en',
        thresh_search   = 2e-6, 
        shared_mem      = 50e9,
        batch_size      = 1,
        matvec          = 4,
        nproc           = None)
    time2 = time.time()

    for rn in range(nroots):
        print("Root:%4d     Var Energy:%12.8f   Gap:%12.8f  CI Dim: %4i "%(rn,e0[rn].real,e0[rn].real-e0[0].real,len(ci_vector)))
    for rn in range(nroots):
        print("Root:%4d     PT  Energy:%12.8f   Gap:%12.8f  CI Dim: %4i "%(rn,e2[rn],e2[rn]-e2[0],len(ci_vector)))

    print("Time spent in the Ex-TPSCI code%16.8f"%(time2-time1))


    for rn in range(nroots):

        e2,_ = compute_pt2_correction(ci_vector[rn], clustered_ham, e0[rn],
                thresh_asci     = 0,
                thresh_search   = 1e-6,
                batch_size      = 10,
                shared_mem          = 2e11,
                pt_type         = 'mp',
                nproc = 24,
                matvec          = 4)

        print("TCI   Root:%4d     Var Energy:%12.8f   Gap:%12.8f  CI Dim: %4i "%(rn,e0[rn].real+ecore,e0[rn].real-e0[0].real,len(ci_vector)))
        print("TCIPT Root:%4d     Var Energy:%12.8f "%(rn,e2+ecore))
