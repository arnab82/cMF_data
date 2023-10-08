import numpy as np
import scipy
import itertools
import time
from math import factorial
import copy as cp
import sys

from pyscf_helper import *
import pyscf
ttt = time.time()
#Hexabenzo[bc,ef,hi,kl,no,qr]coronene
###     PYSCF INPUT
r0 = 2.0
molecule = '''
C   -0.4839   -1.3371    0.0004 
C   -1.3999   -0.2495    0.0004    
C    0.9161   -1.0876    0.0006    
C   -0.9160    1.0876    0.0006    
C    1.3999    0.2495    0.0004    
C    0.4839    1.3372    0.0003    
C   -0.9679   -2.6747   -0.0002    
C   -2.8004   -0.4991    0.0002    
C    1.8324   -2.1756    0.0009    
C   -1.8324    2.1756    0.0009    
C    2.8004    0.4991    0.0002    
C    0.9679    2.6747   -0.0003    
C   -0.0515   -3.7647   -0.0001    
C  -3.7187    0.5893    0.0007   
C  -2.3697   -2.9257   -0.0011    
C  -3.2861   -1.8377   -0.0007  
C   3.2345   -1.9270    0.0010  
C  -1.3489    3.5151    0.0008  
C   1.3490   -3.5151    0.0011  
C  -3.2346    1.9269    0.0012  
C   3.2860    1.8378   -0.0004  
C   2.3697    2.9258   -0.0009   
C   3.7187   -0.5893    0.0003   
C   0.0515    3.7647   -0.0004   
C  -0.5744   -5.0792   -0.0017   
C  -5.1030    0.2975    0.0003   
C  -2.8092   -4.2704   -0.0031   
C  -4.6861   -2.0423   -0.0008   
C   4.1116   -3.0371    0.0017   
C  -2.2938    4.5681    0.0018   
C   2.2938   -4.5681    0.0022   
C  -4.1116    3.0370    0.0020   
C   4.6860    2.0422   -0.0005   
C   2.8092    4.2705   -0.0025   
C   5.1031   -0.2975    0.0000   
C   0.5744    5.0793   -0.0022   
C  -1.9301   -5.3335   -0.0033 
C  -5.5842   -0.9953   -0.0005  
C   3.6540   -4.3385    0.0022  
C   -3.6540    4.3383    0.0023   
C    5.5842    0.9953   -0.0006   
C   1.9302    5.3336   -0.0031    
H   0.0613   -5.9589   -0.0020    
H   -5.8619    1.0735    0.0005    
H   -3.8607   -4.5396   -0.0046   
H   -5.1301   -3.0327   -0.0012    
H    5.1913   -2.9264    0.0019    
H  -2.0011    5.6133    0.0022    
H   2.0011   -5.6133    0.0028    
H  -5.1914    2.9263    0.0025    
H   5.1300    3.0326   -0.0008    
H   3.8607    4.5397   -0.0033    
H   5.8619   -1.0736    0.0001    
H  -0.0614    5.9590   -0.0030   
H  -2.2995   -6.3541   -0.0048    
H  -6.6528   -1.1857   -0.0007    
H   4.3532   -5.1686    0.0029   
H  -4.3533    5.1685    0.0030   
H   6.6528    1.1858   -0.0008  
H   2.2995    6.3543   -0.0045   
''' 
charge = 0
spin  = 0
basis_set = 'sto-3g'

npoly = 21
na = npoly
nb = npoly

### 
orb_basis = 'PM'
cas_nel = 2*npoly
cas_norb = 2*npoly

#Integrals from pyscf
import pyscf
from pyscf import gto, scf, ao2mo, lo
from pyscf.tools import molden
pyscf.lib.num_threads(1)  #with degenerate states and multiple processors there can be issues
#PYSCF inputs

mol = gto.Mole()
mol.atom = molecule

mol.max_memory = 1000 # MB
mol.symmetry = True
mol.charge = charge
mol.spin = spin
mol.basis = basis_set
mol.build()
print("symmertry")
print(mol.topgroup)

#SCF 

#mf = scf.RHF(mol).run(init_guess='atom')
mf = scf.RHF(mol).run()
#C = mf.mo_coeff #MO coeffs
enu = mf.energy_nuc()

if mol.symmetry == True:
    from pyscf import symm
    mo = symm.symmetrize_orb(mol, mf.mo_coeff)
    osym = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mo)
    #symm.addons.symmetrize_space(mol, mo, s=None, check=True, tol=1e-07)
    for i in range(len(osym)):
        print("%4d %8s %16.8f"%(i+1,osym[i],mf.mo_energy[i]))

mo_occ = mf.mo_occ>0
mo_vir = mf.mo_occ==0
print(mo_occ)
print(mo_vir)
mo_occ = np.where(mf.mo_occ>0)[0]
mo_vir = np.where(mf.mo_occ==0)[0]
print(mo_occ)
print(mo_vir)

from pyscf.tools import mo_mapping
s_pop = mo_mapping.mo_comps('C 2pz', mol, mf.mo_coeff)
print(s_pop)
cas_list = s_pop.argsort()[-cas_norb:]
print('cas_list', np.array(cas_list))
print('s population for active space orbitals', s_pop[cas_list])


focc_list = list(set(mo_occ)-set(cas_list))
print(focc_list)
fvir_list = list(set(mo_vir)-set(cas_list))
print(fvir_list)

def get_eff_for_casci(focc_list,cas_list,h,g):
# {{{

    cas_dim = len(cas_list)
    const = 0
    for i in focc_list:
        const += 2 * h[i,i]
        for j in focc_list:
            const += 2 * g[i,i,j,j] -  g[i,j,i,j]

    eff = np.zeros((cas_dim,cas_dim))
    for L,l in enumerate(cas_list):
        for M,m in enumerate(cas_list):
            for j in focc_list:
                eff[L,M] += 2 * g[l,m,j,j] -  g[l,j,j,m]
    return const, eff
# }}}

def reorder_integrals(idx,h,g):
# {{{
    h = h[:,idx] 
    h = h[idx,:] 

    g = g[:,:,:,idx] 
    g = g[:,:,idx,:] 
    g = g[:,idx,:,:] 
    g = g[idx,:,:,:] 
    return h,g
# }}}

local = False
local = True
if local:
    cl_c = mf.mo_coeff[:, focc_list]
    cl_a = lo.Boys(mol, mf.mo_coeff[:, cas_list]).kernel(verbose=4)
    cl_v = mf.mo_coeff[:, fvir_list]
    C = np.column_stack((cl_c, cl_a, cl_v))
else:
    cl_c = mf.mo_coeff[:, focc_list]
    cl_a = mf.mo_coeff[:, cas_list]
    cl_v = mf.mo_coeff[:, fvir_list]
    C = np.column_stack((cl_c, cl_a, cl_v))

from pyscf.tools import mo_mapping
s_pop = mo_mapping.mo_comps('C 2pz', mol, C)
print(s_pop)
cas_list = s_pop.argsort()[-cas_norb:]
print('cas_list', np.array(cas_list))
print('s population for active space orbitals', s_pop[cas_list])
focc_list = list(set(mo_occ)-set(cas_list))
print(focc_list)
fvir_list = list(set(mo_vir)-set(cas_list))
print(fvir_list)



h = C.T.dot(mf.get_hcore()).dot(C)
g = ao2mo.kernel(mol,C,aosym='s4',compact=False).reshape(4*((h.shape[0]),))
const,eff = get_eff_for_casci(focc_list,cas_list,h,g)


focc_list = list(set(mo_occ)-set(cas_list))
print(focc_list)
fvir_list = list(set(mo_vir)-set(cas_list))
print(fvir_list)

ecore = enu + const
h,g = reorder_integrals(cas_list,h,g)
h = h + eff
C = C[:,cas_list]


print("ecore %16.8f"%ecore)


if local:
    print("MULLIKEN")
    m1 = mulliken_ordering(mol,h.shape[0],C)

    ppp = np.column_stack(np.where(m1>.90))
    print(ppp)
    idx = list(ppp[:,1])
    print(m1.shape)
    print(m1)
    m2 = np.where(m1>.90)
    print(m2)
    print(m2[1])

    idx = m2[1]
    import numpy as np
    C = C[:,idx]
    h,g = reorder_integrals(idx,h,g)
    molden.from_mo(mol, 'hbc2.molden', C)
    np.save("ints_h0_hbc2", ecore)
    np.save("ints_h1_hbc2", h)
    np.save("ints_h2_hbc2", g)
    np.save("mo_coeffs_hbc2", C)




