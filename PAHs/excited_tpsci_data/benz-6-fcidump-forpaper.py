import sys, os
import numpy as np
import scipy
import itertools
import time
from math import factorial
import copy as cp
import sys
import tools 

from fermicluster import *
from pyscf_helper import *
import pyscf
from pyscf import gto, scf, ao2mo, molden, lo, mo_mapping, mcscf

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


