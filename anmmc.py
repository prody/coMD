from prody import *
from numpy import *
from random import random
import os.path
import sys
import time

time.sleep(10)
ar = []
for arg in sys.argv:
    ar.append(arg)

initial_pdbn=ar[1]
final_pdbn=ar[2]
initial_pdb_id = initial_pdbn[:initial_pdbn.rfind('.')]
final_pdb_id = final_pdbn[:final_pdbn.rfind('.')]

original_initial_pdb = ar[3]
original_final_pdb = ar[4]

comd_cycle_number = ar[5]

if len(ar) > 6 and ar[6].strip() is not '0':
    devi = float(ar[6])
else:
    devi = 0.5

if len(ar) > 7 and ar[7].strip() is not '0':
    stepcutoff=float(ar[7])
else:
    stepcutoff=2.

if len(ar) > 8 and ar[8].strip() is not '0':
    acceptance_ratio = float(ar[8])
else:
    acceptance_ratio = 0.9

accept_para = 0.1

if len(ar) > 9 and ar[9].strip() is not '0':
    anm_cut=float(ar[9])
else:
    anm_cut=15

if len(ar) > 10 and ar[10].strip() is not '0':
    N=int(ar[10])
else:
    N=10000

if len(ar) > 11 and ar[11].strip() is not '0':
    final_structure_dcd_name = ar[11]
else:
    final_structure_dcd_name = 'cycle_{0}_'.format(int(comd_cycle_number)) + \
                                initial_pdb_id + '_' + final_pdb_id + '_final_structure.dcd'

if len(ar) > 12 and ar[12].strip() is not '0':
    usePseudoatoms = int(ar[12])
else:
    usePseudoatoms = 0

initial_pdb = parsePDB(initial_pdbn)
final_pdb = parsePDB(final_pdbn)

if usePseudoatoms:
    initial_pdb_ca = initial_pdb
    final_pdb_ca = final_pdb
else:
    initial_pdb_ca = initial_pdb.select('name CA or name BB')
    final_pdb_ca = final_pdb.select('name CA or name BB')

# ANM calculation based on current
pdb_anm = ANM('pdb ca')
pdb_anm.buildHessian(initial_pdb_ca, cutoff=anm_cut)
pdb_anm.calcModes()

# Cumulative sum vector preparation for metropolis sampling
eigs = 1/sqrt(pdb_anm.getEigvals())
eigs_n = zeros(eigs.shape)
eigs_n = eigs / sum(eigs)
eigscumsum = eigs_n.cumsum()
U = pdb_anm.getEigvecs()

scale_factor = ones((U.shape[-1]))
if original_initial_pdb != original_final_pdb:
    # Take a step along mode 1 (ID 0) to calculate the scale factor
    pdb_ca = initial_pdb_ca
    pdb_ca_temp = pdb_ca.copy()
    ID = 0
    direction = 1.
    coords_temp = pdb_ca_temp.getCoords()
    coords_temp[0:,0] = coords_temp[0:,0] + direction * U[range(0,len(U),3),ID] * eigs[ID]
    coords_temp[0:,1] = coords_temp[0:,1] + direction * U[range(1,len(U),3),ID] * eigs[ID]
    coords_temp[0:,2] = coords_temp[0:,2] + direction * U[range(2,len(U),3),ID] * eigs[ID]
    pdb_ca_temp.setCoords(coords_temp)
    pdb_ca = pdb_ca_temp.copy()
    biggest_rmsd = calcRMSD(pdb_ca.getCoords(), initial_pdb_ca.getCoords())
    for ID in range(U.shape[-1]):
        scale_factor[ID] = devi/biggest_rmsd # This means that devi is the maximum deviation in RMSD for any step
else: 
    for ID in range(U.shape[-1]):
        pdb_ca = initial_pdb_ca
        pdb_ca_temp = pdb_ca.copy()
        direction = 1.
        coords_temp = pdb_ca_temp.getCoords()
        coords_temp[0:,0] = coords_temp[0:,0] + direction * U[range(0,len(U),3),ID] * eigs[ID] * scale_factor[ID]
        coords_temp[0:,1] = coords_temp[0:,1] + direction * U[range(1,len(U),3),ID] * eigs[ID] * scale_factor[ID]
        coords_temp[0:,2] = coords_temp[0:,2] + direction * U[range(2,len(U),3),ID] * eigs[ID] * scale_factor[ID]
        pdb_ca_temp.setCoords(coords_temp)
        pdb_ca = pdb_ca_temp.copy()
        biggest_rmsd = calcRMSD(pdb_ca.getCoords(), initial_pdb_ca.getCoords())
        scale_factor[ID] = devi/biggest_rmsd # This means that devi is the exact deviation in RMSD for every step

# counts for metropolis sampling
count1 = 0 # Up-hill moves
count2 = 0 # Accepted up-hill moves
count3 = 0 # Down-hill moves

# read MC parameter from file
if os.path.isfile(initial_pdb_id + '_ratio.dat') and os.stat(initial_pdb_id + '_ratio.dat').st_size != 0:
    MCpara = loadtxt(initial_pdb_id + '_ratio.dat')
    accept_para = MCpara[4]
    if MCpara[1] > acceptance_ratio + 0.05:
        accept_para *= 1.5
    elif MCpara[1] < acceptance_ratio - 0.05:
        accept_para /= 1.5
    else:
        savetxt(initial_pdb_id + '_status.dat',[1])
#else:
#    accept_para = 0.1
# MC parameter 1 is the acceptance ratio, which should converge on
# the selected value with a tolerance of 0.05 either side
# and accept_para is adjusted to help bring it within these limits.
# This also happens every 5 steps during the run (lines 173 to 181).

if original_initial_pdb != original_final_pdb:
    # difference from the target structure is defined as the energy and the minimum is zero. 
    native_dist = buildDistMatrix(final_pdb_ca)
    dist = buildDistMatrix(initial_pdb_ca)
    Ep = sum((native_dist - dist)**2)

# Reset pdb_ca (the current structure whole the steps back to the original)
pdb_ca = initial_pdb_ca

step_count = 0
check_step_counts = [0]

sys.stdout.write(' '*2 + 'rmsd' + ' '*2 + 'rand' + ' '*2 + 'ID' + ' '*3 + 'step' \
                 + ' '*2 + 'accept_para' +  ' '*5 + 'f' + '\n')

# MC Loop 
for k in range(N):
    pdb_ca_temp = pdb_ca.copy()
    rand = random()
    ID = argmax(rand<eigscumsum)
    direction = 2*(random()>0.5)-1

    coords_temp = pdb_ca_temp.getCoords()
    coords_temp[0:,0] = coords_temp[0:,0] + direction * U[range(0,len(U),3),ID] * eigs[ID] * scale_factor[ID]
    coords_temp[0:,1] = coords_temp[0:,1] + direction * U[range(1,len(U),3),ID] * eigs[ID] * scale_factor[ID]
    coords_temp[0:,2] = coords_temp[0:,2] + direction * U[range(2,len(U),3),ID] * eigs[ID] * scale_factor[ID]
    pdb_ca_temp.setCoords(coords_temp)

    if original_initial_pdb != original_final_pdb:   
        dist = buildDistMatrix(pdb_ca_temp)
        En = sum((native_dist - dist)**2)

        # Check whether you are heading the right way and accept uphill moves 
        # depending on the Metropolis criterion. Classically this depends on RT 
        # but this is subsumed by the unknown units from having a uniform 
        # spring constant that is set to 1.
        if Ep > En:
            count3 += 1
            pdb_ca = pdb_ca_temp.copy()
            Ep = En
            accepted = 1

        elif exp(-(En-Ep) * accept_para) > random():
            pdb_ca = pdb_ca_temp.copy()
            count1 += 1
            count2 += 1
            Ep = En
            accepted = 1

        else:
            count1 += 1
            accepted = 0

        if count1 == 0:
            f = 1.
        else:
            f = float(count2)/float(count1)

        if (mod(k,5)==0 and not(k==0)):
            # Update of the accept_para to keep the MC para reasonable
            # See comment lines 82 to 85. 
            if f > acceptance_ratio + 0.05:
                accept_para /= 1.5;
            elif f < acceptance_ratio - 0.05:
                accept_para *= 1.5

        if accept_para < 0.001: accept_para = 0.001

    else:
        # for exploration based on one structure
        # all moves are uphill but will be accepted anyway
        pdb_ca = pdb_ca_temp.copy()
        count3 += 1
        accepted = 1
        f = 1.

    rmsd = calcRMSD(pdb_ca.getCoords(), initial_pdb_ca.getCoords())
    sys.stdout.write('{:6.2f}'.format(rmsd) + ' ' + '{:5.2f}'.format(rand) + \
                     '{:4d}'.format(ID) + '{:7d}'.format(k) + ' '*2 + str(accepted) + ' '*2 + \
                     '{:5.4f}'.format(accept_para) + ' '*2 + '{:5.4f}'.format(f) + '\n')

    if rmsd > stepcutoff:
        break
    
# Build an ensemble for writing the final structure to a dcd file
ensemble_final = Ensemble()
ensemble_final.setAtoms(initial_pdb_ca)
ensemble_final.setCoords(initial_pdb_ca)
ensemble_final.addCoordset(pdb_ca.getCoords())
writeDCD(final_structure_dcd_name, ensemble_final)

ratios = [count2/N, count2/count1 if count1 != 0 else 0, count2, k, accept_para ]
savetxt(initial_pdb_id + '_ratio.dat', ratios, fmt='%.2e')

