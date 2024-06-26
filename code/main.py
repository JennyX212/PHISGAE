import numpy as np
import scipy.sparse as sp
import tensorflow as tf
import gc
import random
from clac_metric import cv_model_evaluate
from utils import *
from model import GCNModel
from opt import Optimizer
from pred import cross_validation_experiment

input_path = input("input=")
phage_sim = np.loadtxt(input_path + 'phage_sim.csv', delimiter=',')
host_sim = np.loadtxt(input_path + 'host_I.csv', delimiter=',')
phage_host_matrix = np.loadtxt(input_path + 'phage_host.csv', delimiter=',')
epoch = 4000
emb_dim = 64
lr = 0.01
adjdp = 0.6
dp = 0.2
simw = 1
result = np.zeros((1, 5), float)
average_result = np.zeros((1, 5), float)
circle_time = 1
for i in range(circle_time):
    result += cross_validation_experiment(phage_host_matrix, phage_sim*simw, host_sim*simw, 120, epoch, emb_dim, dp, lr, adjdp)  ###
average_result = result / circle_time
print(average_result)
output_path = input("output=")
np.savetxt(output_path + 'average_result.csv', average_result)
