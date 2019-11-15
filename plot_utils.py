#!/usr/bin/env python
# coding: utf-8

# In[1]:


from statistics import mean
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# In[ ]:





# In[2]:


def n_thresholds(alist, percents=[95], display=True):
    ###Given a list calculates thresholds.
    #   defaults to calculating the 95% threshold
    #   optionally prints the thresholds calculated
    #   returns a dictionary with the thresholds
    #       including and excluding zeros for each

    alist = sorted(alist, reverse=True)
    with_zeros = {}
    for i in percents:
        p = (100.0-float(i))/100.0
        t = float(alist[math.ceil(float(len(alist))*p)])
        with_zeros[i] = t

        if display: print("{0}% threshold: {1}".format(i, t))

    if display: print("\nIgnoring Zeros: ")
    alist = [x for x in alist if (x!=0)]
    skip_zeros = {}
    for i in percents:
        p = (100.0-float(i))/100.0
        t = float(alist[math.ceil(float(len(alist))*p)])
        skip_zeros[i] = t
        if display: print("{0}% threshold: {1}".format(i, t))

    r = {
        'with_zeros':with_zeros,
        'without_zeros':skip_zeros
    }

    return r


# In[3]:


def reproducibility(data):
    num_channels = len(data.columns)
    fig, (axs) = plt.subplots(nrows=num_channels, ncols=num_channels)
    for i_label, i in zip(data.columns.values, range(0,num_channels)):
        for j_label, j in zip(data.columns.values, range(0,num_channels)):
            i_data = data[i_label]
            j_data = data[j_label]
            points = ([],[])#these will represent proteins
            #axs[i][j].set_title("{0} \nvs {1}".format(i_label,j_label))
            for protein in i_data.index:
                 if protein in j_data.index:
                    #here we match proteins from each set
                    points[0].append(i_data[protein])
                    points[1].append(j_data[protein])
            axs[i][j].scatter(points[0],points[1])
    fig.set_figheight(2*num_channels)
    fig.set_figwidth(2*num_channels)
    plt.tight_layout()

