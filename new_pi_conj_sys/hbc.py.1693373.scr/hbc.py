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
#Hexabenzo[a,d,g,j,m,p]coronene
###     PYSCF INPUT
r0 = 2.0
molecule = '''
C    1.4159    0.1095    0.0001    
C   0.6131    1.2808   -0.0001    
C   0.8026   -1.1714   -0.0001    
C   -0.8027    1.1714   -0.0002    
C   -0.6131   -1.2808    0.0002    
C   -1.4158   -0.1095    0.0001    
C    2.8238    0.2182    0.0003    
C    1.2229    2.5546   -0.0001   
C    1.6008   -2.3363   -0.0005   
C   -1.6009    2.3362   -0.0007   
C   -1.2230   -2.5544    0.0003   
C   -2.8237   -0.2183    0.0004   
C    3.5720   -0.9275    0.4330   
C    2.5892    2.6299   -0.4323   
C    3.3874    1.4651   -0.4320    
C   -2.9623    2.2005   -0.4335    
C   -2.5889   -2.6298    0.4336    
C    0.9828   -3.5569   -0.4336    
C   -3.5720    0.9269   -0.4328    
C   -0.4250   -3.6657   -0.4329    
C    2.9624   -2.2011    0.4323    
C   -3.3870   -1.4650    0.4338    
C    0.4246    3.6659    0.4324    
C  -0.9832    3.5571    0.4318    
C    4.8948   -0.8317    0.9518    
C    3.1676    3.8237   -0.9504    
C    4.7096    1.5736   -0.9499    
C   -4.8949    0.8306   -0.9510    
C   -3.1670   -3.8234    0.9526    
C    3.7170   -3.2920    0.9503    
C   -0.9924   -4.8645   -0.9515    
C   -3.7169    3.2910   -0.9528    
C    0.9915    4.8649    0.9511    
C    1.7273   -4.6541   -0.9531    
C   -4.7086   -1.5729    0.9530    
C  -1.7281    4.6548    0.9500    
C    4.4807    3.9039   -1.4076    
C    5.0264   -3.1697    1.4086    
C   -0.2316   -5.9369   -1.4110    
C    5.2584    2.7690   -1.4072    
C    5.6204   -1.9289    1.4094    
C    1.1400   -5.8307   -1.4118    
C   -5.6206    1.9274   -1.4099   
C   -4.4794   -3.9031    1.4116    
C    0.2305    5.9376    1.4090    
C   -5.0265    3.1682   -1.4109    
C   -5.2570   -2.7681    1.4118    
C  -1.1412    5.8318    1.4083    
H    5.3784    0.1359    1.0693   
H   2.5713    4.7263   -1.0677  
H    5.3362    0.6918   -1.0666  
H   -5.3786   -0.1370   -1.0672  
H    3.2664   -4.2756    1.0666  
H   -2.0695   -4.9662   -1.0678  
H    2.0685    4.9665    1.0686  
H   -2.5707   -4.7261    1.0692    
H    2.8071   -4.5888   -1.0704    
H   -3.2662    4.2744   -1.0703    
H   -5.3350   -0.6909    1.0699    
H   -2.8080    4.5898    1.0664    
H    4.8698    4.8403   -1.7946    
H    5.5545   -4.0353    1.7956    
H   -0.7170   -6.8268   -1.7985    
H    1.7562   -6.6354   -1.8001    
H    6.6256   -1.7978    1.7970    
H    6.2721    2.7939   -1.7939    
H   -4.8683   -4.8392    1.7994    
H   -6.6258    1.7959   -1.7971    
H   -5.5545    4.0333   -1.7988    
H    0.7157    6.8277    1.7966    
H   -6.2702   -2.7927    1.7998    
H   -1.7576    6.6367    1.7954   
'''
charge = 0
spin  = 0
basis_set = 'sto-3g'

npoly = 24
na = npoly
nb = npoly

###     TPSCI BASIS INPUT
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

#local = False
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

    ppp = np.column_stack(np.where(m1>.10))
    print(ppp)
    idx = list(ppp[:,1])
    print(m1.shape)
    print(m1)
    m2 = np.where(m1>.10)
    print(m2)
    print(m2[1])

    idx = m2[1]
    import numpy as np
    C = C[:,idx]
    h,g = reorder_integrals(idx,h,g)
    molden.from_mo(mol, 'hbc.molden', C)
    np.save("ints_h0_hbc", ecore)
    np.save("ints_h1_hbc", h)
    np.save("ints_h2_hbc", g)
    np.save("mo_coeffs_hbc", C)



