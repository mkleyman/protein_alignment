from scipy import stats, interpolate, optimize
import math
import numpy as np
from random import shuffle


def curve_distance(ref_times, ref_expression, spline):
    '''
    :param ref_times: time of reference expression
    :param ref_expression: expression values of reference,
     must be in same order as ref_times
    :param spline: spline of expression to compare reference to
    :return: distance between the reference expression and its spline comparison
    '''
    error = 0
    #print "value=" + str(spline(5))
    #print "value=" + str(spline)
    for i, time in enumerate(ref_times):
        #print "spline "+str(np.isnan(spline(time)))
        dif = spline(time) - ref_expression[i]
        error += math.pow(dif, 2)

    return error


def curve_spearman(ref_times, ref_expression, spline):
    '''
       :param ref_times: time of reference expression
       :param ref_expression: expression values of reference,
        must be in same order as ref_times
       :param spline: spline of expression to compare refrence to
       :return: spearman correlation of the refrence expression and
       its spline comaprison
       '''
    spline_expression = [spline(time) for time in ref_times]
    cor, pval = stats.spearmanr(a=ref_expression, b=spline_expression)
    return cor


def curve_pearson(ref_times, ref_expression, spline):
    spline_expression = [spline(time) for time in ref_times]
    cor, pval = stats.pearsonr(ref_expression, spline_expression)
    return cor


def count_above(ref_times, ref_expression, spline, dif_function,min_val):
    if dif_function(ref_times, ref_expression, spline) >= min_val:
        return 1
    else: return 0


def error_function(alignment_function, homolog_dict,
                   ref_dict, dif_function, times,
                   spline_dict, coef):
    '''

    :param alignment_function: function to align reference time series to the
    comparison spline
    :param homolog_dict: dictionary between the reference and comparison
    :param spline_dict: dictionary of the comparison proteins
     and their associated splines
    :param times: time points for reference
    :param ref_dict: dictionary mapping reference protein to expression values
    :param dif_function: function to be used as the difference between 2 time series
    :param coef: arguments for the dif_function
    :return: total error
    '''

    aligned_times = [alignment_function(float(time), coef) for time in times]
    #interval = [aligned_times[0],aligned_times[-1]]
    #spline_dict = make_spline_dict(comp_table, comp_times, interval)
    #ref_proteins = list(ref_table.index)
    tot_error = 0
    for protein in ref_dict:
        tot_error+= dif_function(ref_times=aligned_times,
                                  ref_expression=list(ref_dict[protein]),
                                  spline=spline_dict[homolog_dict[protein]])
    return tot_error


def error_function_summarize(alignment_function, homolog_dict,
                   ref_dict, dif_function, times,
                   spline_dict, summary_fun, weight_dict, coef):
    '''

    :param alignment_function: function to align reference time series to the
    comparison spline
    :param homolog_dict:  dictionary between the reference and comparison list
    :param ref_dict: dictionary mapping reference protein to expression values
    :param dif_function: function to be used as the difference between 2 time series
    :param times:  time points for reference
    :param spline_dict: dictionary of the comparison proteins
     and their associated splines
    :param summary_fun: function of how to summarize error for 1-many mapping of
    homologs between reference an comparison
    :param coef: arguments for the dif_function
    :return: total error
    '''
    aligned_times = [alignment_function(float(time), coef) for time in times]
    # interval = [aligned_times[0],aligned_times[-1]]
    # spline_dict = make_spline_dict(comp_table, comp_times, interval)
    # ref_proteins = list(ref_table.index)
    tot_error = 0
    for protein in ref_dict:
        tot_error += summary_fun([dif_function(ref_times=aligned_times,
                                           ref_expression=list(ref_dict[protein]),
                                           spline=spline_dict[homolog])*weight_dict[homolog] for
                              homolog in homolog_dict[protein]])
    return tot_error


def make_spline_dict(comp_table, times):
    '''

    :param comp_table: pandas table of comparison protein expression
    :param times: time series for comparison
    :return: dictionary of reference proteins and their spline curves
    '''
    proteins = list(comp_table.index)
    times_uniq = list(set(times))
    times_uniq.sort()
    #str_times = list(comp_table.values)
    #str_times.pop()
    #times = [float(time) for time in str_times]
    spline_dict = {}
    weight_dict = {}
    for protein in proteins:
        #start = min(interval[0], times[0])
        #end = max(interval[1], times[-1])
        vals,weight = median_values(times, list(comp_table.loc[protein]))
        spline_dict[protein] = interpolate.UnivariateSpline(x=times_uniq,
                                     y=vals, s=50)
                                     #bbox=[0, end], check_finite=True)
        weight_dict[protein] = weight
        if np.isnan(spline_dict[protein](5)):
            print "problem please average your data"
        #print "val at 5:"+str(spline_dict[protein](5.0))
        #print "val at 3000:" + str(spline_dict[protein](3000))
        #print times
        #print "y-vals: "+str(list(comp_table.loc[protein]))


    return spline_dict,weight_dict

def make_spline_dict_rand(comp_table, times):
    '''

    :param comp_table: pandas table of comparison protein expression
    :param times: time series for comparison
    :return: dictionary of reference proteins and their spline curves
    '''
    proteins = list(comp_table.index)
    times_uniq = list(set(times))
    times_uniq.sort()
    spline_dict = {}
    for protein in proteins:
        #start = min(interval[0], times[0])
        #end = max(interval[1], times[-1])
        vals = median_values(times, list(comp_table.loc[protein]))
        shuffle(vals)
        spline_dict[protein] = interpolate.UnivariateSpline(x=times_uniq,
                                     y=vals)
                                     #bbox=[0, end], check_finite=True)
        if np.isnan(spline_dict[protein](5)):
            print "problem please average your data"


    return spline_dict



def mean_values(times, values):
    '''

    :param times: times in time series, this will include repeats
    :param values: values for each time for a single protein
    :return: mean values for each time point
    '''
    exp_dict = {}
    for i,time in enumerate(times):
        if time in exp_dict:
            exp_dict[time].append(values[i])
        else:  exp_dict[time] = [values[i]]
    pairs = exp_dict.items()
    pairs.sort(key=lambda tup: tup[0])
    return [float(sum(pair[1]))/float(len(pair[1])) for pair in pairs]


def median_values(times, values):
    '''
    :param times: times in time series, this will include repeats
    :param values: values for each time for a single protein
    :return: median values for each time point
    '''
    exp_dict = {}
    for i,time in enumerate(times):
        if time in exp_dict:
            exp_dict[time].append(values[i])
        else:  exp_dict[time] = [values[i]]
    pairs = exp_dict.items()
    pairs.sort(key=lambda tup: tup[0])
    within_time_var = np.mean([np.var(pair[1]) for pair in pairs])
    values = [np.median(pair[1]) for pair in pairs]
    return values, np.var(values)/within_time_var


def make_ref_dict(ref_table,times, avg_function):
    '''

    :param ref_table: pandas table references
    :param times: times in refrences time series
    :param avg_function: function to summarize between values at same time point
    :return: dictionary of proteins to "average" expression at each time point
    '''
    ref_dict= {}
    ref_proteins = list(ref_table.index)
    for protein in ref_proteins:
        ref_dict[protein] = avg_function(times,list(ref_table.loc[protein]))
    return ref_dict


