#!/usr/bin/env python
# coding: utf-8

from statistics import mean, stdev, variance
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
    

def skipZeroMean(l):
    l = [x for x in l if (x!=0)]
    if len(l)==0: return 0
    return mean(l)

def cal_fold_changes_and_variances(samples, approxZero=None):
    #Calculates the fold changes and variances for the graphs.
    #  Called by the graph_by_variance function.
    #   Optionally, specify a number to replace an average of zero with.
    #   Using zero to calculate fold change gives error or zero,
    #   yet a full on-off change is noteworthy and should be marked.
    #   By default, uses have the minimum non-zero value
    variances = {}
    fold_changes = {}
    
    if approxZero == None: localZero = True
    else: localZero = False
        
    sample_names = list(samples.keys())
    
    for ser_index in range(0,len(sample_names)-1): #these are the keys to samples
        ser = sample_names[ser_index]
        sample_df1 = samples[ser]
        for o_ser_index in range(ser_index+1,len(sample_names)): #keys again
            variances_averaged = {} #used in graph to show v1 vs v2
            o_ser = sample_names[o_ser_index]
            #compare variance in sample versus otherSample
            sample_df2 = samples[o_ser]
            if localZero: 
                approxZero = get_least_value({ ser:sample_df1,o_ser:sample_df2 }) / 2.0
            for protein in sample_df1.index:
                t1 = [x for x in sample_df1.loc[protein,:] if x != 0]
                t2 = [x for x in sample_df2.loc[protein,:] if x != 0]

                if len(t1) > 1: v1 = variance(t1)
                else: v1 = 0
                if len(t2) > 1: v2 = variance(t2)
                else: v2 = 0
                v = (v1+v2)/2
                variances_averaged[v1]=v2
                variances[((ser,o_ser),protein)] = v
                #This implements an approximate zero to avoid divided by zero errors
                m2 =skipZeroMean(t2); m1 = skipZeroMean(t1)
                if m1 ==0: m1 = approxZero
                if m2 ==0: m2 = approxZero
                fold_changes[((ser,o_ser),protein)] = m2/m1

    returnValues = {
        "fold_changes":fold_changes,
        "variances":variances
    }
        
    return returnValues   
    
def graph_by_variance(data, cell_types = ["Inter","Motor"], 
                      FOLD_CHANGE_THRESHOLD = 2, approxZero = None,
                      threshold95 = None, threshold99 = None):
    #Graphs the variance vs fold change plot.
    #  Takes the data and replicate dict.
    #  By default, marks 2 fold doubling/halving lines
    #  By default, uses half the minimum non-zero in 
    #   calculating fold change in place of zero 
    #   to avoid divide by zero errors.
    #  Marks the 95% variance threshold.
    #  Optionally marks the 99% threshold.
    
    #data = data.apply(np.log)
    data.fillna(0, inplace=True)
    
    #Splits into samples
    samples = {}
    r_names=np.array(data.columns.values)
    for cell_type in cell_types:
        cells_of_type = list(s for i,s in enumerate(r_names) if cell_type in s)
        data_by_type = data.loc[:, cells_of_type]
        samples[cell_type] = data_by_type
    
    if approxZero == None:            
        mins = []
        for l in samples:
            l_mins = samples[l].apply(
                lambda l_list: min([x for x in l_list if (x!=0)])).tolist()
            for i in l_mins: mins.append(i)
        approxZero = min(mins) / 2.0
        
    #Calculates points
    fc_v = cal_fold_changes_and_variances(samples, approxZero)
    fold_changes= fc_v["fold_changes"]
    variances = fc_v["variances"]
    
    if threshold95 ==None:
        variances_sorted = sorted(variances.values(), reverse=True)
        threshold95 = float(variances_sorted[math.ceil(float(len(variances_sorted))*.05)])
        #print (threshold95)
    
    #Volcano Graph of all the data
    plt.rc('axes', titlesize=55)
    plt.rc('axes', labelsize=30)
    plt.rc('xtick', labelsize=25)
    plt.rc('ytick', labelsize=25)
    
    fig= plt.figure(figsize=(24,16))
    
    plt.title("Log 2 Fold Changes")
    
    log2_fold_changes = [math.log2(x) for x in fold_changes.values()]
    plt.scatter(log2_fold_changes, variances.values())
    plt.axvline(x=math.log2(FOLD_CHANGE_THRESHOLD), linestyle='dashed')
    plt.axvline(x=-math.log2(FOLD_CHANGE_THRESHOLD), linestyle='dashed')
    plt.axhline(y=threshold95, color='b', linestyle='-')
    if threshold99 != None:
        plt.axhline(y=threshold99, color='b', linestyle='-',alpha=.5)
    plt.gca().invert_yaxis()
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("Technical Variance")
    plt.show()

def get_altered(data, cell_types = ["Inter","Motor"],
                FOLD_CHANGE_THRESHOLD = 2, VAR_THRES = 95, approxZero=None):
    
    #data = data.apply(np.log)
    data.fillna(0, inplace=True)
    
    #Splits into samples
    samples = {}
    r_names=np.array(data.columns.values)
    for cell_type in cell_types:
        cells_of_type = list(s for i,s in enumerate(r_names) if cell_type in s)
        data_by_type = data.loc[:, cells_of_type]
        samples[cell_type] = data_by_type
        
    if approxZero == None:            
        mins = []
        for l in samples:
            l_mins = samples[l].apply(
                lambda l_list: min([x for x in l_list if (x!=0)])).tolist()
            for i in l_mins: mins.append(i)
        approxZero = min(mins) / 2.0
        
    fc_v = cal_fold_changes_and_variances(samples, approxZero)
    fold_changes= fc_v["fold_changes"]
    fd = pd.DataFrame.from_dict(fold_changes, orient='index')
    fd = fd.sort_values(0)
    variances = fc_v["variances"]
    
    variances_sorted = sorted(variances.values(), reverse=True)
    variances_for_thresh = [x for x in variances_sorted if x !=0]
    threshold = float(variances_sorted[math.ceil(float(len(variances_for_thresh))*((100.0-VAR_THRES)/100.0))])
    
    sig_proteins = []
    for key in fd.index:
        if abs(math.log2(fold_changes[key])) > math.log2(FOLD_CHANGE_THRESHOLD):
            if variances[key] < threshold:
                sig_proteins.append(key[1])
                                       
    df = data.loc[sig_proteins, :]
    return df