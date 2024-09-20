import numpy as np
import sys
import os

# force params
m = 50
k = 10
ol = 0.25
d = 1

#############################
# main parameters to modify #
#############################

fld = sys.argv[1]
sim = int(sys.argv[2])
# mcParams = np.loadtxt(fld+'/params.csv', delimiter=',')
# params = mcParams[sim, :]
pST = int(sys.argv[3]) #phenotype state transition 
dp = float(sys.argv[4])
kp = float(sys.argv[5])

print("starting simulation with phenotype state transition: {}, death probability factor: {}, kill probability factor: {}".format(pST, dp, kp) )

cancerPDL1_m = 0.01 
tcellMigBias = 0.076
cd4Diff = 0.132
cd4PDL1 = 0.139
macMigBias = 0.123
macM1 = 0.187
macM2 = 0.433
macPDL1 = 0.264 
cd8RecRate = 0 #0.001 #0.0001 #ensure just cancer cells
cd4RecRate = 0 #0.001 #ensure just cancer cells
macRecRate = 0 # 0.006 #ensure just cancer cells
recDelay = 2.85
necroticGrowth = 0 #params[12]
necrosisLimit = 940.3

p_ferroptosis_sensitive = 0.5 #probability that a cell is ferroptosis sensitive

#############################
# ------------------------- #
#############################

cellParams = np.zeros((13, 4))

# cancer params
cellParams[0, 0] = m  # mu
cellParams[1, 0] = k  # kc
cellParams[2, 0] = d  # damping
cellParams[3, 0] = ol  # overlap
cellParams[4, 0] = 20.0  # diameter (um)
cellParams[5, 0] = 0.02 #1/24.0  # div probability (hours) # Gong 2017
cellParams[6, 0] = 0 #1/(24.0*5.0) # death probability (hours) # Gong 2017
cellParams[7, 0] = 40.0  # influence distance
cellParams[8, 0] = cancerPDL1_m  # pdl1 when expressed
cellParams[9, 0] = 1e6  # prob of gaining pdl1 (is multiplied by t cell influence) -> pretty sure the cancer cells from the images already had pdl1 when injected (check with evanthia). 1e6 pretty much guarantees that cancer cells will express pdl1
cellParams[10, 0] = p_ferroptosis_sensitive #probability a cell is ferroptotis sensitive 
cellParams[11, 0] = 50 #ferroptosis influence distance


# cd4 params
cellParams[0, 1] = m  # mu
cellParams[1, 1] = k  # kc
cellParams[2, 1] = d  # damping
cellParams[3, 1] = ol  # overlap
cellParams[4, 1] = 10.0  # diameter (um)
cellParams[5, 1] = 1/(24.0*3.0) # lifespan (days) # Gong 2017 (for CD8 cells)
cellParams[6, 1] = 90.0  # migration speed base um/hr
cellParams[7, 1] = cd4Diff  # differentiation to Treg
cellParams[8, 1] = 40.0  # influence distance
cellParams[9, 1] = tcellMigBias  # migration bias base
cellParams[10, 1] = cd4PDL1  # pdl1 (when Treg)

# cd8 params
cellParams[0, 2] = m  # mu
cellParams[1, 2] = k  # kc
cellParams[2, 2] = d  # damping
cellParams[3, 2] = ol  # overlap
cellParams[4, 2] = 10.0  # diameter (um)
cellParams[5, 2] = 1/(24.0*3.0*dp) # lifespan (days) Gong 2017
cellParams[6, 2] = 90.0  # migration speed um/hr | https://onlinelibrary.wiley.com/doi/epdf/10.1038/icb.2012.75 -> scaled based on model scale
cellParams[7, 2] = 0.1*kp # killProb Gong 2017 (other works use very different probabilities
cellParams[8, 2] = 2.0  # infScale -> arbitrarily set
cellParams[9, 2] = 40.0  # influence distance
cellParams[10, 2] = tcellMigBias  # migration bias base
cellParams[11, 2] = 0.053  # proliferation prob

"""
this parameterizes the rate at which a modeled T cell will progress through phenotypic state transitions
as articulated by the embedded Gene Regulatory Network. A higher value results in a faster progression through states. 
This value must be an integer. A recommended range is 1-10 but should not exceed the length of the smallest trajectory
"""
cellParams[12, 2] = pST

# macrophage params
cellParams[0, 3] = m  # mu
cellParams[1, 3] = k  # kc
cellParams[2, 3] = d  # damping
cellParams[3, 3] = ol  # overlap
cellParams[4, 3] = 20.0  # diameter (um)
cellParams[5, 3] = 1/(24.0*3.0) # (set based on cd8) lifespan (days) -> https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4040600/ | Petty and Yang, Tumor-associated macrophages: implications in cancer immunotherapy, 2017
cellParams[6, 3] = 20.0  # migration speed um/hr | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702045/ -> scaled based on model scale
cellParams[7, 3] = macM1  # kM1
cellParams[8, 3] = macM2  # kM2
cellParams[9, 3] = 40.0  # influence distance
cellParams[10, 3] = macMigBias  # migration bias
cellParams[11, 3] = macPDL1  # pdl1

recParams = np.zeros((5, 1))
recParams[0] = cd8RecRate # cd8RecRate
recParams[1] = macRecRate # mRecRate
recParams[2] = cd4RecRate # cd4RecRate
recParams[3] = 50.0  # recDist (recruit a uniform distribution recDist away from the tumor edge
#recParams[3] = 0.25 # max recruitment cytokine conc
recParams[4] = recDelay # recruitment delay (days)

envParams = np.zeros((7, 1))
envParams[0] = 5.0  # initTumorSize x | circle radius
envParams[1] = 18 # simulation duration (days)
envParams[2] = necroticGrowth # necrotic growth
envParams[3] = 0.5 # necrotic region outward force
envParams[4] = necrosisLimit # necrosis limit (accounts for diffusion limit of oxygen, but is adjustable based on the scale of the simulation)
envParams[5] =  15 #ferroptosis delay (days)
envParams[6] = 50 #ferroptosis radius (um) 

# os.system('mkdir -p ' + sys.argv[1] + '/params')
saveFld = sys.argv[1]+'/set_'+sys.argv[2]+'/params'
os.system('mkdir -p '+saveFld)
np.savetxt(saveFld+'/cellParams.csv', cellParams, delimiter=',')
np.savetxt(saveFld+'/recParams.csv', recParams, delimiter=',')
np.savetxt(saveFld+'/envParams.csv', envParams, delimiter=',')

# np.savetxt(sys.argv[1] + '/params/cellParams.csv', cellParams, delimiter=',')
# np.savetxt(sys.argv[1] + '/params/recParams.csv', recParams, delimiter=',')
# np.savetxt(sys.argv[1] + '/params/envParams.csv', envParams, delimiter=',')


print('------Parameters Grid Done-----\n')
