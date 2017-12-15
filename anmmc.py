from prody import *
from numpy import *
from random import random
import os.path
import sys

ar =[]
for arg in sys.argv:
        ar.append(arg)

initial_pdbn=ar[1]
final_pdbn=ar[2]

if len(ar) > 3:
    anm_cut=float(ar[3])
else:
    anm_cut=15

if len(ar) > 4:
    devi=float(ar[4])
else:
    devi=0.1

if len(ar) > 5:
    accept_para=int(ar[5])
else:
    accept_para=0.1

if len(ar) > 6:
    N=int(ar[6])
else:
    N=1000000

original_initial_pdb = ar[-2]
original_final_pdb = ar[-1]

initial_pdb_id = initial_pdbn.split('.')[0]
final_pdb_id = final_pdbn.split('.')[0]

fo = open('anmmc_log.txt','w')

initial_pdb = parsePDB(initial_pdbn)
final_pdb = parsePDB(final_pdbn)

# Current Structure 
initial_pdb_ca = initial_pdb.ca
stepcutoff = 0.5 * (len(initial_pdb_ca) ** 0.5)

# ANM calculation based on current
if os.path.isfile(initial_pdb_id + '.anm.npz'):
    pdb_anm = loadModel(initial_pdb_id + '.anm.npz')
    fo.write("ANM model loaded\n")
else:
    # ANM calculation based on current
    pdb_anm = ANM('pdb ca')
    # Build Hessian Matrix
    pdb_anm.buildHessian(initial_pdb_ca, cutoff=anm_cut)
    pdb_anm.calcModes(n_modes='all')
    saveModel(pdb_anm,initial_pdb_id,matrices=True)

# Cumulative sum vector preparation for metropolis sampling
eigs = 1/sqrt(pdb_anm.getEigvals())
eigs_nn = zeros(eigs.shape)
eigs_nn[:] = eigs[:]
eigs = eigs / sum(eigs)
eigscumsum = eigs.cumsum()
U = pdb_anm.getEigvecs()

# Target Structure
final_pdb_ca = final_pdb.ca

# Number of residues on protein structure         
size=initial_pdb_ca.getResnums().shape[0]    

# Difference between current and final structure 
deviation = final_pdb_ca.getCoords() - initial_pdb_ca.getCoords()

# Scale factor for steps
scale_factor = sqrt(abs(devi*min(pdb_anm.getEigvals())))

# counts for metropolis sampling
count1 = 0 # Up-hill moves
count2 = 0 # Accepted up-hill moves
count3 = 0 # Down-hill moves

# read MC parameter from file
if os.path.isfile(initial_pdb_id + '_ratio.dat') and os.stat(initial_pdb_id + '_ratio.dat').st_size != 0:
    MCpara = loadtxt(initial_pdb_id + '_ratio.dat')
    accept_para = MCpara[4]
    if MCpara[1] > 0.95:
        accept_para *= 1.5
    elif MCpara[1] < 0.85:
        accept_para /= 1.5
    else:
        savetxt(initial_pdb_id + '_status.dat',[1])
else:
    accept_para = 0.1
# The best value for MC parameter 1 is around 0.9 
# so values less than 0.85 and greater than 0.95 are not preferred 
# and accept_para is adjusted to help bring it within these limits.
# This also happens every 25 steps during the run (lines 144 to 150).
    
# difference from the target structure is defined as the energy and the minimum is zero. 
native_dist = buildDistMatrix(final_pdb_ca)
dist = buildDistMatrix(initial_pdb_ca)
Ep = sum((native_dist - dist)**2)

pdb_ca = initial_pdb_ca

ensemble = Ensemble()
ensemble.setAtoms(initial_pdb_ca)
ensemble.setCoords(initial_pdb_ca)

ensemble_final = Ensemble()
ensemble_final.setAtoms(initial_pdb_ca)
ensemble_final.setCoords(initial_pdb_ca)

step_count = 0
check_step_counts = [0]
#exit
# MC Loop 
for k in range(N):
    pdb_ca_temp = pdb_ca.copy()
    rand = random()    
    ID = argmax(rand<eigscumsum)        
    direction = 2*(random()>0.5)-1

    coords_temp = pdb_ca_temp.getCoords()
    coords_temp[0:,0] = coords_temp[0:,0] + direction * U[range(0,len(U),3),ID] * eigs_nn[ID] * scale_factor
    coords_temp[0:,1] = coords_temp[0:,1] + direction * U[range(1,len(U),3),ID] * eigs_nn[ID] * scale_factor
    coords_temp[0:,2] = coords_temp[0:,2] + direction * U[range(2,len(U),3),ID] * eigs_nn[ID] * scale_factor
    pdb_ca_temp.setCoords(coords_temp)
   
    dist = buildDistMatrix(pdb_ca_temp)
    En = sum((native_dist - dist)**2)

    if original_initial_pdb != original_final_pdb:
        # check whether you are heading the right way
        # and accept uphill moves depending on the
        # Metropolis criterion
        # The energies depend on the scale_factor,
        # which acts like an effective temperature.
        if Ep > En:
            count3 += 1
            pdb_ca = pdb_ca_temp.copy()
            Ep = En
            step_count += 1
        elif exp(-(En-Ep)*accept_para) > random():
            pdb_ca = pdb_ca_temp.copy() 
            count1 += 1
            count2 += 1
            Ep = En
            step_count += 1
        else:
            count1 += 1

        if (mod(k,25)==0 and not(k==0)):
            # Update of the accept_para to keep the MC para reasonable
            # See comment lines 82 to 85. 
            if count2/count1 > 0.95:
                accept_para*=1.5;
            elif count2/count1 < 0.85:
                accept_para/=1.5

    else:
        # for exploration based on one structure (two runs)
        # all moves are uphill but will be accepted anyway
        pdb_ca = pdb_ca_temp.copy()
        count1 += 1
        count2 += 1
        Ep = En
        step_count += 1

    if (mod(k,1000)==0 and not(k==0)):
        check_step_counts.append(step_count)
        writeDCD(initial_pdb_id + '_' + final_pdb_id + '_steps_' + str(check_step_counts[-2]) + '-' \
                 + str(check_step_counts[-1]-1) + '.dcd',ensemble[check_step_counts[-2]:check_step_counts[-1]-1])

    coord_diff = pdb_ca.getCoords() - initial_pdb_ca.getCoords()
    fo.write(str(En) + '\t' + str(linalg.norm(coord_diff.ravel())) + '\t' + str(rand) + '\t' + str(ID) + '\t' + str(k) + '\t' + str(step_count) '\n')

    if linalg.norm(coord_diff.ravel()) > stepcutoff: 
        break
        
    ensemble.addCoordset(pdb_ca.getCoords())
    
ensemble_final.addCoordset(pdb_ca.getCoords())
    
writeDCD(initial_pdb_id + '_' + final_pdb_id + '_final_structure.dcd', ensemble_final)
writeDCD(initial_pdb_id + '_' + final_pdb_id + '_ensemble.dcd', ensemble)
ratios = [count2/N, count2/count1 if count1 != 0 else 0, count2, k, accept_para ]
savetxt(initial_pdb_id + '_ratio.dat', ratios)

