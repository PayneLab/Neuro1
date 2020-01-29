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

def reproducibility(data, pair=None, logging=True, show_only=[]):
    #Plots points from a dataframe against two columns
    #    for example, the proteins in sample 1 against sample 2
    #
    #  data: a pandas dataframe
    #  pair: a list of two (or more) to graph;
    #    defaults to None, which graphs all columns
    #  logging: whether the data should be graphed on a log-scale
    #  show_only: a list of points to show; a subset from data
    #  color_subset=(color, list)

    if pair:
        num_channels = len(pair)
        channels = pair
    else:
        num_channels = len(data.columns.values)
        channels = data.columns.values

    if show_only != []:
        data=data.loc[show_only]

    if num_channels == 2: num_graph_rows = 1
    else: num_graph_rows = num_channels

    plt.rcParams['figure.dpi'] = 200
    fig, (axs) = plt.subplots(nrows=num_graph_rows, ncols=num_graph_rows)
    scs = []
    annots = []

    for i_label, i in zip(channels, range(0,num_graph_rows)):
        scs.append([])
        annots.append([])
        i_data = data[i_label]
        if num_graph_rows == 1: channels = channels[i+1:]
        for j_label, j in zip(channels, range(0,num_graph_rows)):
            j_data = data[j_label]
            points = ([],[],[])#these will represent proteins

            for protein in i_data.index:
                 if protein in j_data.index:
                    #here we match proteins from each set
                    points[0].append(j_data[protein])
                    points[1].append(i_data[protein])
                    points[2].append(protein)

            if num_graph_rows > 1: #If we are graphing multiple
                ax = axs[i][j]
                #only hide ticks if we are graphing multiple
                ax.set_xticks(ticks=[])#hides ticks
                ax.set_yticks(ticks=[])#hides ticks
                #If these are the edge graphs
                if j == 0: ax.set_ylabel(i_label, rotation=45,horizontalalignment='right')
                if i == num_graph_rows-1: ax.set_xlabel(j_label, rotation=45,horizontalalignment='right')
            else:
                #if only two are being compared, the labels are oriented normally
                ax = axs #psuedonym for consistency with multiple graphs
                ax.set_ylabel(i_label)
                ax.set_xlabel(j_label)

            sc = ax.scatter(points[0],points[1])

            #forces the scales equal
            plt.gca().set_aspect('equal', adjustable='box')
            limit = max(ax.get_ylim()[1],ax.get_xlim()[1]) *1.001# add a small buffer
            if logging: min_limit = 90
            #elif not logging: min_limit = -1
            ax.set_ylim(min_limit, limit)
            ax.set_xlim(min_limit, limit)
            scs[i].append(sc)

            annot = ax.annotate("", xy=(0,0), xytext=(20,20),
                    textcoords="offset points", bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"))
            annot.set_visible(False)
            annots[i].append(annot)

            if logging:
                ax.set_yscale('log')
                ax.set_xscale('log')

            if num_graph_rows > 1:
                #only hide ticks if we are graphing multiple
                ax.set_xticks(ticks=[])#hides ticks
                ax.set_yticks(ticks=[])#hides ticks

    def update_annot(ind, annot, sc):
        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = "{}".format("/n".join([points[2][n] for n in ind["ind"]]))
        annot.set_text(text)
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(event):
        for i in range(0,num_graph_rows):
            for j in range(0,num_graph_rows):
                if num_graph_rows > 1:
                    ax = axs[i][j]
                else: ax = axs
                if event.inaxes == ax:
                    cont, ind = scs[i][j].contains(event)
                    if cont:
                        annots[i][j].set_visible(True)
                        update_annot(ind, annots[i][j], scs[i][j])
                        fig.canvas.draw_idle()
                    else:
                        if annots[i][j].get_visible():
                            annots[i][j].set_visible(False)
                            fig.canvas.draw_idle()
    fig.canvas.mpl_connect("motion_notify_event", hover)

    fig.set_figheight(2*num_channels)
    fig.set_figwidth(2*num_channels)
    plt.show()
