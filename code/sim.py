import numpy as np
import pandas as pd

from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics.pairwise import pairwise_distances


############# #############   cosine_similarity   ############# #############
phage_f = np.loadtxt('../data/dnapro_f/phage_dna_pro.csv', delimiter=',')
phage_sim =  cosine_similarity(phage_f)
print(phage_sim.shape)

np.savetxt("../data/dnapro_f/phage_dnapro_sim.csv", phage_sim, delimiter=",")

############# #############   GIP_similarity   ############# #############
import math

def Gussian_similarity(intMat):

    gamall = 1

    nl = np.shape(intMat)[0]

    print(nl)

    # calculate gamal for Gaussian kernel calculation
    sl = np.zeros(nl)

    for i in range(nl):
        sl[i] = np.square(np.linalg.norm(intMat[i, :]))
    gamal = nl / sum(np.transpose(sl)) * gamall
    print(gamal)

    # hostMat = np.zeros([nl,nl],float)
    # for i in range(nl):
    #     for j in range(nl):
    #         hostMat[i, j] = math.exp(-gamal*np.square(np.linalg.norm(intMat[:, i]-intMat[:, j])))
    # return hostMat

    phageMat = np.zeros([nl, nl], float)
    for i in range(nl):
        for j in range(nl):
            phageMat[i, j] = math.exp(-gamal * np.square(np.linalg.norm(intMat[i, :] - intMat[j, :])))

    return phageMat

if __name__ == "__main__":
    intMat = np.loadtxt("../data/dnapro_f/phage_dna_pro.csv", delimiter=',')
    # hostMat = Gussian_similarity(intMat)
    phageMat = Gussian_similarity(intMat)

np.savetxt('../data/dnapro_f/phage_gipsim.csv', phageMat, delimiter=",")





