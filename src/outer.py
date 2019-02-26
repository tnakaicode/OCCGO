import numpy as np
import matplotlib.pyplot as plt


def v_inner(v, w):
    num = v.shape
    val = 0
    for i in range(num[0]):
        val += v[i] * w[i]
    return val


def v_normal(v):
    num = v.shape
    val = v_inner(v, v)
    return v / np.sqrt(val)


def v_outer(v, w):
    num = v.shape
    index = num[0]
    array = []
    for i in range(index):
        iv, iw = (i+1) % index, (i+2) % index
        #print (iv, iw)
        val = v[iv]*w[iw] - v[iw]*w[iv]
        array.append(val)
    func = np.array(array)
    return func
