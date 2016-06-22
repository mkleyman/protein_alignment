import parser2
import align
import curves
from scipy import optimize
import numpy as np
import csv, random
import significance, opt_params


reference = "../processed/mouse_proteins.csv"
comparison = "../processed/human_proteins.csv"
homolog_file = "../raw_data/HOM_MouseHumanSequence.rpt"
#go_file = "../raw_data/MOUSE_10090_idmapping_selected.tab"
#plot_folder = "../plots"

mouse_proteins, times = parser2.read_in_expression(reference)
#print mouse_proteins.columns.values
#print set(mouse_proteins.index)
print mouse_proteins.shape
#print times

#print mouse_proteins[:5]
#print list(mouse_proteins.iloc[1])


human_proteins, htimes = parser2.read_in_expression(comparison)
print human_proteins.shape
#print human_proteins.columns.values

ref_dict, comp_dict = parser2.read_in_homologs(filename=homolog_file, homolog_col=0,
                                               organism_col=1, protein_col=12,
                                               reference="mouse, laboratory",
                                               comparison="human")
print len(ref_dict.keys())
print ref_dict.items()[1]
print len(comp_dict.keys())
print comp_dict.items()[1]

homolog_dict,homologs = parser2.reference_to_comparison_listmap(ref_dict, comp_dict,
                                                                list(mouse_proteins.index),
                                                                set(human_proteins.index))
print len(homolog_dict.keys())
print homolog_dict.items()[1]
#print list(human_proteins.loc['Q04917'])
spline_dict,weight_dict = align.make_spline_dict(human_proteins, htimes)
'''for key in spline_dict:
    plotting.plot_spline(htimes, spline_dict[key],plot_folder,key)'''
#normalize weight dict
avg = np.mean(weight_dict.values())
for key in weight_dict:
    weight_dict[key] = 1.0
    #weight_dict[key] = weight_dict[key]/avg

print len(spline_dict.keys())
#rint spline_dict.items()[1]

mouse_subset = mouse_proteins[mouse_proteins.index.isin(homolog_dict)]
human_subset = human_proteins[human_proteins.index.isin(homologs)]
mouse_dict = align.make_ref_dict(mouse_subset, times, align.mean_values)
times_uniq = list(set(times))
times_uniq.sort()
#parser.read_in_GO(go_file, set(homolog_dict.keys()))
negative_spearman = lambda ref_times, ref_expression, spline: \
    -1.0*align.curve_spearman(ref_times, ref_expression, spline)
negative_pearson = lambda ref_times, ref_expression, spline: \
    -1.0 *align.curve_pearson(ref_times, ref_expression, spline)
print "key# "+str(len(mouse_dict.keys()))
count_spearman = lambda ref_times, ref_expression, spline: \
    align.count_above(ref_times,ref_expression,spline,align.curve_spearman, 0.6)
count_pearson =  lambda ref_times, ref_expression, spline: \
    align.count_above(ref_times,ref_expression,spline,align.curve_pearson, 0.8)
#print "init_err: "+str(align.error_function_summarize(curves.linear, homolog_dict, mouse_dict,
 #                          count_pearson, times_uniq, spline_dict, max, weight_dict, [1.1*100,.01]))

error_fun = lambda a: len(mouse_dict.keys())-align.error_function_summarize(curves.linear,
                                           homolog_dict, mouse_dict,
                                           count_spearman, times_uniq,
                                            spline_dict,max, weight_dict,a)
error_fun2 = lambda a: align.error_function_summarize(curves.linear,
                                           homolog_dict, mouse_dict,
                                            align.curve_distance, times_uniq, spline_dict,min,a)
error_fun3 = lambda a: len(mouse_dict.keys()) - align.error_function_summarize(curves.exp,
                                                    homolog_dict, mouse_dict,
                                                    count_spearman, times_uniq, spline_dict, max, a)
error_fun4 = lambda a: len(mouse_dict.keys()) - align.error_function_summarize(curves.logit,
                                                                               homolog_dict, mouse_dict,
                                                                               count_spearman, times_uniq, spline_dict,
                                                                               max, a)
error_fun5 = lambda a: len(mouse_dict.keys()) - align.error_function_summarize(curves.linear,
                                                                             homolog_dict, mouse_dict,
                                                                           count_pearson, times_uniq, spline_dict,
                                                                             max, weight_dict,a)
error_fun6 = lambda a: align.error_function_summarize(curves.linear,
                                           homolog_dict, mouse_dict,
                                           count_spearman, times_uniq,
                                            spline_dict,max, weight_dict,a)
error_fun7 = lambda a: align.error_function_summarize(curves.linear,
                                                      homolog_dict, mouse_dict,
                                                      count_pearson, times_uniq,
                                                      spline_dict, max,
                                                      weight_dict, a)
error_fun8 = lambda a: align.error_function_summarize(curves.exp2,
                                                      homolog_dict, mouse_dict,
                                                      count_pearson, times_uniq,
                                                      spline_dict, max,
                                                      weight_dict, a)
error_fun9 = lambda a: align.error_function_summarize(curves.linear,
                                                      homolog_dict, mouse_dict,
                                                      align.curve_spearman, times_uniq,
                                                      spline_dict, max,
                                                      weight_dict, a)

error_fun10 = lambda a: align.error_function_summarize(curves.linear,
                                                      homolog_dict, mouse_dict,
                                                      align.curve_pearson, times_uniq,
                                                      spline_dict, max,
                                                      weight_dict, a)
error_fun11 = lambda a: align.error_function_summarize(curves.linear,
                                                      homolog_dict, mouse_dict,
                                                      align.curve_distance, times_uniq,
                                                      spline_dict, max,
                                                      weight_dict, a)
error_funq = lambda a: len(mouse_dict.keys()) - align.error_function_summarize(curves.quadratic,
                                                                              homolog_dict, mouse_dict,
                                                                              count_pearson, times_uniq,
                                                                              spline_dict, max, weight_dict, a)

error_fun_pearson_linear = lambda a: align.error_function_summarize(curves.linear,
                                      homolog_dict, mouse_dict,
                                      negative_pearson, times_uniq,
                                      spline_dict, max, weight_dict, a)
error_fun_pearson_quadratic = lambda a: align.error_function_summarize(curves.quadratic,
                                            homolog_dict, mouse_dict,
                                            negative_pearson, times_uniq,
                                            spline_dict, max, weight_dict, a)

error_fun_pearson_cubic = lambda a: align.error_function_summarize(curves.quadratic,
                                       homolog_dict, mouse_dict,
                                       negative_pearson, times_uniq,
                                       spline_dict, max, weight_dict, a)
error_fun_pearson_count_cubic = lambda a: -1.0*align.error_function_summarize(curves.cubic2,
                                                       homolog_dict, mouse_dict,
                                                       count_pearson, times_uniq,
                                                       spline_dict, max, weight_dict, a)

error_fun_pearson_count_sqrt = lambda a: -1.0 * align.error_function_summarize(curves.square_root,
                                                                        homolog_dict, mouse_dict,
                                                                        count_pearson, times_uniq,
                                                                        spline_dict, max, weight_dict, a)

error_fun_pearson_count_log = lambda a: -1.0 * align.error_function_summarize(curves.logarithm,
                                                                               homolog_dict, mouse_dict,
                                                                               count_pearson, times_uniq,
                                                                               spline_dict, max, weight_dict, a)
#result = optimize.basinhopping(error_fun, [1,0.05], niter=300, minimizer_kwargs=options)
#result = optimize.minimize(error_fun, [1,0.05], method="Nelder-Mead")
'''
rranges = (slice(0,3300,100),slice(1,300,10))

range1 = xrange(0,3300,100)
range2 = xrange(0,180,5)
recording = param_search.record_twoparam(error_fun10, range1,range2)
np.savetxt("linear_param_pearson.txt", recording)
#recording = np.loadtxt("linear_param_search.txt")
heatmap = sns.heatmap(recording, cmap="YlGnBu",xticklabels=10, yticklabels= 10)
heatmap.get_figure().savefig("linear_param_pearson.png")
'''


'''
result = optimize.brute(error_fun,rranges, full_output=True, finish= optimize.fmin)
#print result.x
print "coefs: "+ str(result[0])
print "best: "+ str(result[1])
#err = align.error_function(curves.linear, homolog_dict, mouse_dict,
                           #align.curve_spearman, times_uniq, spline_dict, [2.09575833e+02,   2.95969535e-02])
count_good = align.error_function_summarize(curves.linear, homolog_dict, mouse_dict,
                           align.curve_spearman, times_uniq, spline_dict, max, weight_dict, result[0])       '''

options = {"method":"Nelder-Mead"}
use_fun = error_fun_pearson_count_log
init_conditions = [300,20,10,50]
brute = False
'''
print "initial err: "+str(use_fun(init_conditions))
if(not brute):
    result = optimize.basinhopping(use_fun,init_conditions, minimizer_kwargs=options)
    print result.x
    print use_fun(result.x)
else:
    rranges = (slice(0.0,3300.0,50.0),slice(1.0,100.0,2.0), slice(0,500,5))
    result = optimize.brute(use_fun,rranges, full_output=True, finish= optimize.fmin)
    print "coefs: " + str(result[0])
    print "best: " + str(result[1])
'''

thresholds = [0.6]
init_dict = {}
functions = {}
functions["Linear"] = curves.linear
init_dict["Linear"] = [2000.0,25.0]
functions["Quadratic"] =curves.quadratic
init_dict["Quadratic"] = [1000.0,25.0,5.0]
#functions["Cubic"] = curves.cubic
#init_dict["Cubic"] = [50.0,25.0,7.0,3.0]
functions["Sqrt"] = curves.square_root
init_dict["Sqrt"] = [2000.0,25.0]

functions["Gompertz"] = curves.gom
init_dict["Gompertz"] = [40.0,5.0,3.5]
options = {"method":"Nelder-Mead"}
funct_list = ["Sqrt","Linear","Gompertz","Quadratic"]
with open("results_partial.csv",'wb') as resultfile:
    result_writer = csv.writer(resultfile)
    random_seeds = [random.random() for i in xrange(50)]
    for seed in random_seeds:
        for functname in funct_list:
            funct = functions[functname]
            init = init_dict[functname]
            for thresh in thresholds:
                print functname
                count_spearman = lambda ref_times, ref_expression, spline: \
                align.count_above(ref_times,ref_expression,spline,align.curve_spearman, thresh)
                error_fun_spearman_count_pars = lambda a: -1.0 * align.error_function_summarize_rand(funct,
                                               homolog_dict, mouse_dict,
                                               count_spearman, times_uniq,
                                               spline_dict, max, weight_dict, seed,a)
                result = optimize.basinhopping(error_fun_spearman_count_pars,init, minimizer_kwargs=options)
                filename = functname+"_spearman_revolve.txt"
                with open(filename,mode='a') as writefile:
                    writefile.write(str(error_fun_spearman_count_pars(result.x)))






                count_pearson =  lambda ref_times, ref_expression, spline: \
                align.count_above(ref_times,ref_expression,spline,align.curve_pearson, thresh)
                error_fun_pearson_count_pars = lambda a: -1.0 * align.error_function_summarize_rand(funct,
                                               homolog_dict, mouse_dict,
                                               count_pearson, times_uniq,
                                               spline_dict, max, weight_dict, seed,a)
                result = optimize.basinhopping(error_fun_pearson_count_pars,init, minimizer_kwargs=options)
                filename = functname+"_pearson_revolve.txt"
                with open(filename,mode='a') as writefile:
                    writefile.write(str(error_fun_pearson_count_pars(result.x)))
