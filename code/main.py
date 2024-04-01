import numpy as np
import scipy.sparse as sp
import tensorflow as tf
import gc
import random
from clac_metric import cv_model_evaluate
from utils import *
from model import GCNModel
from opt import Optimizer


def PredictScore(train_phage_host_matrix, phage_matrix, host_matrix, seed, epochs, emb_dim, dp, lr,  adjdp):
    np.random.seed(seed)
    tf.reset_default_graph()
    tf.set_random_seed(seed)
    adj = constructHNet(train_phage_host_matrix, phage_matrix, host_matrix)
    adj = sp.csr_matrix(adj)
    association_nam = train_phage_host_matrix.sum()
    X = constructNet(train_phage_host_matrix)

    # phage_kmer = np.loadtxt('../data/phage_kmer_order.csv', delimiter=',')
    # host_kmer = np.loadtxt('../data/host_kmer_order.csv', delimiter=',')
    # X = np.vstack((phage_kmer, host_kmer))

    features = sparse_to_tuple(sp.csr_matrix(X))
    num_features = features[2][1]
    features_nonzero = features[1].shape[0]
    adj_orig = train_phage_host_matrix.copy()
    adj_orig = sparse_to_tuple(sp.csr_matrix(adj_orig))

    adj_norm = preprocess_graph(adj)
    adj_nonzero = adj_norm[1].shape[0]
    placeholders = {
        'features': tf.sparse_placeholder(tf.float32),
        'adj': tf.sparse_placeholder(tf.float32),
        'adj_orig': tf.sparse_placeholder(tf.float32),
        'dropout': tf.placeholder_with_default(0., shape=()),
        'adjdp': tf.placeholder_with_default(0., shape=())
    }
    model = GCNModel(placeholders, num_features, emb_dim,
                     features_nonzero, adj_nonzero, train_phage_host_matrix.shape[0], name='PHISGAE')
    with tf.name_scope('optimizer'):
        opt = Optimizer(
            preds=model.reconstructions,
            labels=tf.reshape(tf.sparse_tensor_to_dense(
                placeholders['adj_orig'], validate_indices=False), [-1]),
            model=model,
            lr=lr, num_u=train_phage_host_matrix.shape[0], num_v=train_phage_host_matrix.shape[1], association_nam=association_nam)
    sess = tf.Session()
    sess.run(tf.global_variables_initializer())

    for epoch in range(epochs):
        feed_dict = dict()
        feed_dict.update({placeholders['features']: features})
        feed_dict.update({placeholders['adj']: adj_norm})
        feed_dict.update({placeholders['adj_orig']: adj_orig})
        feed_dict.update({placeholders['dropout']: dp})
        feed_dict.update({placeholders['adjdp']: adjdp})
        _, avg_cost = sess.run([opt.opt_op, opt.cost], feed_dict=feed_dict)
        if epoch % 100 == 0:
            feed_dict.update({placeholders['dropout']: 0})
            feed_dict.update({placeholders['adjdp']: 0})
            res = sess.run(model.reconstructions, feed_dict=feed_dict)
            print("Epoch:", '%04d' % (epoch + 1),
                  "train_loss=", "{:.5f}".format(avg_cost))
    print('Optimization Finished!')
    feed_dict.update({placeholders['dropout']: 0})
    feed_dict.update({placeholders['adjdp']: 0})
    res = sess.run(model.reconstructions, feed_dict=feed_dict)
    sess.close()
    return res

# def updateX(train_phage_host_matrix, phage_matrix, host_matrix, seed, epochs, emb_dim, dp, lr,  adjdp):
#     np.random.seed(seed)
#     tf.reset_default_graph()
#     tf.set_random_seed(seed)
#     adj = constructHNet(train_phage_host_matrix, phage_matrix, host_matrix)
#     adj = sp.csr_matrix(adj)
#     association_nam = train_phage_host_matrix.sum()
#     X = constructHNet(train_phage_host_matrix, phage_matrix, host_matrix)
#     features = sparse_to_tuple(sp.csr_matrix(X))
#     num_features = features[2][1]
#     features_nonzero = features[1].shape[0]
#     adj_orig = train_phage_host_matrix.copy()
#     adj_orig = sparse_to_tuple(sp.csr_matrix(adj_orig))
#
#     adj_norm = preprocess_graph(adj)
#     adj_nonzero = adj_norm[1].shape[0]
#     placeholders = {
#         'features': tf.sparse_placeholder(tf.float32),
#         'adj': tf.sparse_placeholder(tf.float32),
#         'adj_orig': tf.sparse_placeholder(tf.float32),
#         'dropout': tf.placeholder_with_default(0., shape=()),
#         'adjdp': tf.placeholder_with_default(0., shape=())
#     }
#     model = GCNModel(placeholders, num_features, emb_dim,
#                      features_nonzero, adj_nonzero, train_phage_host_matrix.shape[0], name='PHISGAE')
#
#     with tf.name_scope('optimizer'):
#         opt = Optimizer(
#             preds=model.reconstructions,
#             labels=tf.reshape(tf.sparse_tensor_to_dense(
#                 placeholders['adj_orig'], validate_indices=False), [-1]),
#             model=model,
#             lr=lr, num_u=train_phage_host_matrix.shape[0], num_v=train_phage_host_matrix.shape[1], association_nam=association_nam)
#     sess = tf.Session()
#     sess.run(tf.global_variables_initializer())
#
#     for epoch in range(epochs):
#         feed_dict = dict()
#         feed_dict.update({placeholders['features']: features})
#         feed_dict.update({placeholders['adj']: adj_norm})
#         feed_dict.update({placeholders['adj_orig']: adj_orig})
#         feed_dict.update({placeholders['dropout']: dp})
#         feed_dict.update({placeholders['adjdp']: adjdp})
#         _, avg_cost = sess.run([opt.opt_op, opt.cost], feed_dict=feed_dict)
#         if epoch % 100 == 0:
#             feed_dict.update({placeholders['dropout']: 0})
#             feed_dict.update({placeholders['adjdp']: 0})
#             res = sess.run(model.reconstructions, feed_dict=feed_dict)
#             print("Epoch:", '%04d' % (epoch + 1),
#                   "train_loss=", "{:.5f}".format(avg_cost))
#     print('Optimization Finished!')
#     feed_dict.update({placeholders['dropout']: 0})
#     feed_dict.update({placeholders['adjdp']: 0})
#     XL_embeddings = sess.run(model.embeddings, feed_dict=feed_dict)
#     sess.close()
#     return XL_embeddings

def cross_validation_experiment(phage_host_matrix, phage_matrix, host_matrix, seed, epochs, emb_dim, dp, lr, adjdp):
    index_matrix = np.mat(np.where(phage_host_matrix == 1))
    association_nam = index_matrix.shape[1]
    random_index = index_matrix.T.tolist()
    random.seed(seed)
    random.shuffle(random_index)
    k_folds = 5
    CV_size = int(association_nam / k_folds)
    temp = np.array(random_index[:association_nam - association_nam %
                                 k_folds]).reshape(k_folds, CV_size,  -1).tolist()
    temp[k_folds - 1] = temp[k_folds - 1] + \
        random_index[association_nam - association_nam % k_folds:]
    random_index = temp
    metric = np.zeros((1, 7))
    print("seed=%d, evaluating phage-host...." % (seed))
    for k in range(k_folds):
        print("------this is %dth cross validation------" % (k+1))
        train_matrix = np.matrix(phage_host_matrix, copy=True)
        train_matrix[tuple(np.array(random_index[k]).T)] = 0
        phage_len = phage_host_matrix.shape[0]
        host_len = phage_host_matrix.shape[1]
        phage_host_res = PredictScore(train_matrix, phage_matrix, host_matrix, seed, epochs, emb_dim, dp, lr,  adjdp)
        predict_y_proba = phage_host_res.reshape(phage_len, host_len)
        # predict_y_proba2 = (predict_y_proba == predict_y_proba.max(axis=1)[:, None]).astype(int)

        # XL_embeddings = updateX( train_matrix, phage_matrix, host_matrix, seed, epochs, emb_dim, dp, lr, adjdp)  ###
        # XLname = "../result/XL" + str(k) + ".csv"
        # np.savetxt(XLname, XL_embeddings, fmt='%f', delimiter=",")
        # predA = "../result/predA" + str(k) + ".csv"
        # np.savetxt(predA, predict_y_proba, fmt='%f', delimiter=",")

        metric_tmp = cv_model_evaluate(
            phage_host_matrix, predict_y_proba, train_matrix)
        print(metric_tmp)
        metric += metric_tmp
        del train_matrix
        gc.collect()
    print(metric / k_folds)
    metric = np.array(metric / k_folds)
    return metric


if __name__ == "__main__":
    phage_sim = np.loadtxt('../data/phage_sim.csv', delimiter=',')
    host_sim = np.loadtxt('../data/host_I.csv', delimiter=',')
    phage_host_matrix = np.loadtxt('../data/phage_host.csv', delimiter=',')
    # PHX = np.loadtxt('../data/X000.csv', delimiter=',')  ######
    epoch = 4000
    emb_dim = 64
    lr = 0.01
    adjdp = 0.6
    dp = 0.2
    simw = 1
    result = np.zeros((1, 7), float)
    average_result = np.zeros((1, 7), float)
    circle_time = 1
    for i in range(circle_time):
        result += cross_validation_experiment(phage_host_matrix, phage_sim*simw, host_sim*simw, 120, epoch, emb_dim, dp, lr, adjdp)  ###
    average_result = result / circle_time
    print(average_result)
