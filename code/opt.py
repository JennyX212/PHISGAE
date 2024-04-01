from clr import cyclic_learning_rate
import tensorflow as tf


class Optimizer():
    def __init__(self, model, preds, labels, lr, num_u, num_v, association_nam):
        norm = num_u*num_v / float((num_u*num_v-association_nam) * 2)
        preds_sub = preds
        labels_sub = labels
        pos_weight = float(num_u*num_v-association_nam)/(association_nam)

        global_step = tf.Variable(0, trainable=False)
        self.optimizer = tf.train.AdamOptimizer(learning_rate=cyclic_learning_rate(global_step=global_step, learning_rate=lr*0.1,
                                                                                   max_lr=lr, mode='exp_range', gamma=.995))
        preds_sub2 = tf.reshape(preds_sub, [num_u, num_v])
        preds_sub = tf.reshape(preds_sub2, [-1])
        cost2 = 1 / (num_u * num_v) * tf.reduce_sum(tf.reduce_max(preds_sub2, axis=1)-tf.reduce_sum(preds_sub2, axis=1))

        rho = 0.002
        lam = 0.02
        l1_regularizer = tf.keras.regularizers.l1(lam)
        penalty = l1_regularizer(preds_sub2)
        penalty = 1 / (num_u * num_v) * penalty

        self.cost = norm * tf.reduce_mean(
            tf.nn.weighted_cross_entropy_with_logits(
                logits=preds_sub, targets=labels_sub, pos_weight=pos_weight))-rho* norm * cost2+ norm * penalty  #


        self.opt_op = self.optimizer.minimize(
            self.cost, global_step=global_step,)
        self.grads_vars = self.optimizer.compute_gradients(self.cost)
