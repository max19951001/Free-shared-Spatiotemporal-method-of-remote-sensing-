# -*- coding: utf-8 -*-

import numpy as np


def initialize_parameters(parameters=None):

    parameters = {} if not parameters else parameters

    def safe_pull_value(parameters, key, default):
        return parameters.get(key, default)

    # number of clusters desired
    K = safe_pull_value(parameters, 'k', 6)

    # maximum number of iterations
    I = safe_pull_value(parameters, 'I', 100)

    # maximum of number of pairs of clusters which can be merged
    P = safe_pull_value(parameters, 'P', 2)

    # threshold value for  minimum number of samples in each cluster(discarding clusters)
    minS = safe_pull_value(parameters, 'THETA_M', 10)

    # threshold value for standard deviation (for split)
    maxStdv = safe_pull_value(parameters, 'maxStdv', 0.1)
    # threshold value for pairwise distances (for merge)
    minDis = safe_pull_value(parameters, 'minDis', 2)

    # percentage of change in clusters between each iteration
    # (to stop algorithm)
    M = 0.05

    # number of starting clusters
    k = safe_pull_value(parameters, 'k', K)

    ret = locals()
    ret.pop('safe_pull_value')
    ret.pop('parameters')
    globals().update(ret)


def myISODATA(img, parameters=None):
    """
    Classify a numpy 'img' using Isodata algorithm.
    Parameters: a dictionary with the following keys.
            - img: an input numpy array that contains the image to classify.
            - parameters: a dictionary with the initial values.
              If 'parameters' are not specified, the algorithm uses the default
              ones.
                  + number of clusters desired.
                    k = 6
                  + max number of iterations.
                    I = 100
                  + max number of pairs of clusters which can be merged.
                    P = 2
                  + threshold value for min number in each cluster.
                    minS = 10
                  + threshold value for standard deviation (for split).
                    maxStdv= 0.1
                  + threshold value for pairwise distances (for merge).
                    minDis = 2
                  + threshold change in the clusters between each iter.
                    M = 0.05
        Note: if some(or all) parameters are not provided, default values
              will be used.
    Returns:
            - img_class: a numpy array with the classification.
    """
    global K, I, P, maxStdv, minDis, minS, M, k
    initialize_parameters(parameters)
    shape = img.shape
    img = img.reshape(-1, shape[2])
    m = np.zeros((shape[0] * shape[1], 1))
    pad = np.zeros(img.shape)
    c = np.zeros((k, shape[2]))

    for i in range(k):
        c[i] = img[i, :]
    for iter in range(I):

        distance = np.zeros((shape[0] * shape[1], k))

        c_sum = np.zeros((k, shape[2]))
        c_num = np.zeros((k, 1))
        cnt = 0
        for i in range(k):
            pad = np.tile(c[i], (shape[0] * shape[1], 1))
            pad = img - pad
            distance[:, i] = np.linalg.norm(pad, axis=1)

        m = np.argmin(distance, axis=1)
        for i in range(k):
            t_result = m - i
            t_value = img[np.argwhere(t_result == 0)]
            c_num[i] = t_value.shape[0]
            c_sum[i] = t_value.sum(axis=0)

        if (c_num < minS).any():
            deleteDic = np.zeros(k)
            t = k
            for i in range(t):
                if c_num[i] < minS:
                    k = k - 1
                    deleteDic[i] = 1
            t_distance = np.zeros((shape[0] * shape[1], k))
            t_c = np.zeros((k, shape[2]))
            t_i = 0
            for i in range(t):
                if deleteDic[i] == 0:
                    t_distance[:, t_i] = distance[:, i]
                    t_c[t_i, :] = c[i]
                    t_i = t_i + 1
            distance = t_distance
            c = t_c
            c_sum = np.zeros((k, shape[2]))
            c_num = np.zeros((k, 1))
            cnt = 0
            m = np.argmin(distance, axis=1)
            for i in range(k):
                t_result = m - i
                t_value = img[np.argwhere(t_result == 0)]
                c_num[i] = t_value.shape[0]
                c_sum[i] = t_value.sum(axis=0)

        if ((iter % 2) == 0) or k < (K / 2):
            b_split = False
            t_maxStd = -1
            for i in range(k):
                t_result = m - i
                t_value = img[np.argwhere(t_result == 0)]
                std = np.std(t_value, axis=0)
                if (std > maxStdv).any():
                    t_n_feature = np.argmax(std)
                    if std[0, t_n_feature] > t_maxStd:
                        t_maxStd = std[0, t_n_feature]
                        n_feature = t_n_feature.copy()
                        n_class = i
                        b_split = True

            if b_split:

                split_t_result = m - n_class
                split_t_value = img[np.argwhere(split_t_result == 0)]
                std = np.std(split_t_value, axis=0)

                k = k + 1
                t_row1 = c[n_class, :]
                t_row2 = t_row1.copy()
                t_row1[n_feature] = t_row1[n_feature] - M * std[0, n_feature]
                t_row2[n_feature] = t_row2[n_feature] + M * std[0, n_feature]

                c[n_class, :] = t_row1.T

                c = np.r_['0,2,1', c, t_row2.T]
                distance = np.zeros((shape[0] * shape[1], k))
                c_sum = np.zeros((k, shape[2]))
                c_num = np.zeros((k, 1))
                cnt = 0
                for i in range(k):
                    pad = np.tile(c[i], (shape[0] * shape[1], 1))
                    pad = img - pad
                    distance[:, i] = np.linalg.norm(pad, axis=1)

                m = np.argmin(distance, axis=1)
                for i in range(k):
                    t_result = m - i
                    t_value = img[np.argwhere(t_result == 0)]
                    c_num[i] = t_value.shape[0]
                    c_sum[i] = t_value.sum(axis=0)
        if ((iter % 2) == 1) or k > (K * 2):
            b_merge = False
            EucildDistence = np.zeros((k, k))
            for classi in range(k):
                t_class = np.tile(c[classi, :], (k, 1))
                t_minus = c - t_class
                EucildDistence[:, classi] = np.linalg.norm(t_minus, axis=1)
            t_height = k * (k - 1) / 2
            t_height = np.uint32(t_height)
            distStruct = np.zeros((t_height, 5))
            cursor = 0
            for classi in range(1, k):
                for classj in range(0, classi):
                    distStruct[cursor, :] = [EucildDistence[classi, classj], classi, classj, 0, 0]
                    cursor = cursor + 1
            distStruct = distStruct[np.lexsort(distStruct[:, ::-1].T)]
            for i in range(t_height):
                if distStruct[i, 4] == 0 and distStruct[i, 0] < minDis:
                    b_merge = True
                    distStruct[i, 3] = 1
                    for j in range(t_height):
                        if distStruct[j, 1] == distStruct[i, 1] or distStruct[j, 2] == distStruct[i, 1] or distStruct[j, 1] == distStruct[i, 2] or distStruct[j, 2] == distStruct[i, 2]:
                            distStruct[j, 4] = 1

            t_c = c.copy()
            marker = False
            for i in range(t_height):
                if distStruct[i, 3] == 1:
                    class_a = distStruct[i, 1]
                    class_b = distStruct[i, 2]
                    class_a = np.uint32(class_a)
                    class_b = np.uint32(class_b)
                    k = k - 1
                    #
                    t_c[class_a, :] = (t_c[class_a, :] + t_c[class_b, :]) / 2
                    t_c[class_b, :] = np.zeros((1, shape[2]))
                    marker = True
            if marker:
                c = t_c[np.nonzero(t_c)]
                c = c.reshape(k, shape[2])
                distance = np.zeros((shape[0] * shape[1], k))
                c_sum = np.zeros((k, shape[2]))
                c_num = np.zeros((k, 1))
                cnt = 0
                for i in range(k):
                    pad = np.tile(c[i], (shape[0] * shape[1], 1))
                    pad = img - pad
                    distance[:, i] = np.linalg.norm(pad, axis=1)

                m = np.argmin(distance, axis=1)
                for i in range(k):
                    t_result = m - i
                    t_value = img[np.argwhere(t_result == 0)]
                    c_num[i] = t_value.shape[0]
                    c_sum[i] = t_value.sum(axis=0)

        m = m.flatten()
        m = np.uint8(m)

        if (c_num == 0).any():
            zero_index = np.argwhere(c_sum == 0)
            c_sum[zero_index] = 0.01
            c_num[zero_index] = 1
        t = c_sum / c_num
        if (t == c).all():
            m = m.reshape((shape[0], shape[1]))
            return m, c
        c = t

    m = m.reshape((shape[0], shape[1]))
    return m, c
