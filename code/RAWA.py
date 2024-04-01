import numpy as np
import pandas as pd

phagenamelist = pd.read_csv('../data/species/p_name.csv', header=None)
hostnamelist = pd.read_csv('../data/species/h_name_s.csv', header=None)
rawA = pd.read_csv('../data/species/rawA_s.csv', header=None)

print(phagenamelist.shape)
print(hostnamelist.shape)
print(rawA.shape)

num_p = phagenamelist.shape[0]
num_h = hostnamelist.shape[0]

newA = np.zeros((num_p,num_h))

num_i = rawA.shape[0]
for i in range(num_i):
    ip = np.array(phagenamelist)[:, 1].tolist().index(rawA.iloc[i, 1])
    ih = np.array(hostnamelist)[:, 1].tolist().index(rawA.iloc[i, 2])
    newA[ip, ih] = 1

# print(newA)
# print(newA.sum())
# print(newA.sum(axis=0))
np.savetxt("../data/species/phage_host.csv", newA, fmt = '%d', delimiter=",")

