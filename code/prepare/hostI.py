
import numpy as np


HH = np.loadtxt("../data/species/h_name_s.csv", delimiter=',')
num_h = np.shape(HH)[0]

hostsim = np.eye(num_h)
np.savetxt("../data/species/host_I.csv", hostsim, fmt = '%d', delimiter=",")

'''
num_h = 293
hostsim = np.zeros((num_h,num_h))
np.savetxt("../data/species/host_sim0.csv", hostsim, fmt = '%d', delimiter=",")
'''