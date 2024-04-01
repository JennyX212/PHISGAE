import numpy as np
import tensorflow as tf
import scipy.sparse as sp


def weight_variable_glorot(input_dim, output_dim, name=""):
    init_range = np.sqrt(6.0/(input_dim + output_dim))
    initial = tf.random_uniform(
        [input_dim, output_dim],
        minval=-init_range,
        maxval=init_range,
        dtype=tf.float32
    )

    return tf.Variable(initial, name=name)


def dropout_sparse(x, keep_prob, num_nonzero_elems):
    noise_shape = [num_nonzero_elems]
    random_tensor = keep_prob
    random_tensor += tf.random_uniform(noise_shape)
    dropout_mask = tf.cast(tf.floor(random_tensor), dtype=tf.bool)
    pre_out = tf.sparse_retain(x, dropout_mask)
    return pre_out*(1./keep_prob)

## 这个函数的作用就是 返回一个稀疏矩阵的非0值坐标、非0值和整个矩阵的shape
def sparse_to_tuple(sparse_mx):
    if not sp.isspmatrix_coo(sparse_mx):
        sparse_mx = sparse_mx.tocoo()
    coords = np.vstack((sparse_mx.row, sparse_mx.col)).transpose()
    values = sparse_mx.data
    shape = sparse_mx.shape
    return coords, values, shape


def preprocess_graph(adj):
    adj_ = sp.coo_matrix(adj)
    rowsum = np.array(adj_.sum(1))
    degree_mat_inv_sqrt = sp.diags(np.power(rowsum, -0.5).flatten())
    adj_nomalized = adj_.dot(
        degree_mat_inv_sqrt).transpose().dot(degree_mat_inv_sqrt)
    adj_nomalized = adj_nomalized.tocoo()
    return sparse_to_tuple(adj_nomalized)


def constructNet(phage_host_matrix):
    phage_matrix = np.matrix(
        np.zeros((phage_host_matrix.shape[0], phage_host_matrix.shape[0]), dtype=np.int8))
    host_matrix = np.matrix(
        np.zeros((phage_host_matrix.shape[1], phage_host_matrix.shape[1]), dtype=np.int8))

    mat1 = np.hstack((phage_matrix, phage_host_matrix))
    mat2 = np.hstack((phage_host_matrix.T, host_matrix))
    adj = np.vstack((mat1, mat2))
    return adj


def constructHNet(phage_host_matrix, phage_matrix, host_matrix):
    mat1 = np.hstack((phage_matrix, phage_host_matrix))
    mat2 = np.hstack((phage_host_matrix.T, host_matrix))
    return np.vstack((mat1, mat2))
