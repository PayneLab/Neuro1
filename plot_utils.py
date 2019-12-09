#!/usr/bin/env python
# coding: utf-8

from statistics import mean
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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

def reproducibility(data, pair=None):
    if pair: 
        num_channels = len(pair)
        channels = pair
    else: 
        num_channels = len(data.columns.values)
        channels = data.columns.values
        
    if num_channels == 2: num_channels = 1
        
    fig, (axs) = plt.subplots(nrows=num_channels, ncols=num_channels)
    scs = []
    annots = []
    
    for i_label, i in zip(channels, range(0,num_channels)):
        scs.append([])
        annots.append([])
        i_data = data[i_label]
        if num_channels == 1: channels = channels[i+1:]
        for j_label, j in zip(channels, range(0,num_channels)):
            j_data = data[j_label]
            points = ([],[],[])#these will represent proteins
            #axs[i][j].set_title("{0} \nvs {1}".format(i_label,j_label))
               
            for protein in i_data.index:
                 if protein in j_data.index:
                    #here we match proteins from each set
                    points[0].append(j_data[protein])
                    points[1].append(i_data[protein])
                    points[2].append(protein)
            if num_channels > 1: #If we are graphing multiple
                ax = axs[i][j]
                ax.axis('equal')
                #only hide ticks if we are graphing multiple
                ax.set_xticks(ticks=[])#hides ticks
                ax.set_yticks(ticks=[])#hides ticks
                #If these are the edge graphs
                if j == 0: ax.set_ylabel(i_label, rotation=45,horizontalalignment='right')
                if i == num_channels-1: ax.set_xlabel(j_label, rotation=45)
            else: ax = axs
            sc = ax.scatter(points[0],points[1])
            scs[i].append(sc)
            
            annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points", bbox=dict(boxstyle="round", fc="w"), arrowprops=dict(arrowstyle="->"))
            annot.set_visible(False)
            annots[i].append(annot)
            
    def update_annot(ind, annot, sc):
        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = "{}".format("/n".join([points[2][n] for n in ind["ind"]]))
        annot.set_text(text)
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(event):
        for i in range(0,num_channels):
            for j in range(0,num_channels):
                if num_channels > 1:
                    ax = axs[i][j]
                else: ax = axs
                if event.inaxes == ax:
                    cont, ind = scs[i][j].contains(event)
                    if cont:
                        update_annot(ind, annots[i][j], scs[i][j])
                        annots[i][j].set_visible(True)
                        fig.canvas.draw_idle()
                    else:
                        if annots[i][j].get_visible():
                            annots[i][j].set_visible(False)
                            fig.canvas.draw_idle()
    fig.canvas.mpl_connect("motion_notify_event", hover)
    
    fig.set_figheight(2*num_channels)
    fig.set_figwidth(2*num_channels)
    #plt.tight_layout()
    plt.show()

