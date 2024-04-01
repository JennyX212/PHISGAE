import os, sys
sys.path.append(os.getcwd())
import random
random.seed(1)
import numpy as np
import pandas as pd
import torch
import torch.autograd as autograd
import torch.nn as nn
import torch.optim as optim
torch.manual_seed(1)
from sklearn.model_selection import LeaveOneOut
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score
import warnings
warnings.filterwarnings("ignore")


def obtainfeatures(data,file_path1,strs):
    host_features=[]
    for i in data:
        host_features.append(np.loadtxt(file_path1+i+strs).tolist())
    return np.array(host_features)

data=pd.read_csv('../data/species/p_name.txt',header=None,sep=',')
#data=data[data[2]==1]

datap = data[0].values.tolist()


#datap = list(set(datap))
#p_name = pd.DataFrame(datap)

#p_name.to_csv('../data/data2/p_name.csv', header = None, encoding='utf-8-sig')

phage_feature_pro=obtainfeatures(datap,'../data/Protein_f/phage_pro_feature/','.txt')
#phage_all=np.concatenate((phage_feature_dna, phage_feature_pro),axis=1)

np.savetxt("../data/Protein_f/phage_pro_feature.csv", phage_feature_pro, delimiter=",")

print(phage_feature_pro.shape)

