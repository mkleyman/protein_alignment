import csv
import ast
import re
import curves
import parser2
import align
import numpy as np
import significance
from scipy import stats
'''
results = []
with open("results.csv",mode='r') as resultfile:
    resultreader = csv.reader(resultfile)
    for row in resultreader:
        result = row
        result[3] = ast.literal_eval(re.sub('[0-9.]\s+',',',result[3]))
        print len(result[3])
        results.append(result)
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
with open("rate.csv",mode='wb') as ratefile:
    ratewriter = csv.writer(ratefile)
    for row in results:
        for mouse_age in xrange(16,49):
            nrow = row[:]
            nrow[4]= functions[row[0]](mouse_age, row[3])/365.0
            nrow[3]= mouse_age
            ratewriter.writerow(nrow)
'''
reference = "../processed/mouse_proteins.csv"
comparison = "../processed/human_proteins.csv"
homolog_file = "../raw_data/HOM_MouseHumanSequence.rpt"
mouse_go_file = "../raw_data/MOUSE_10090_idmapping_selected.tab"
human_go_file = "../raw_data/HUMAN_9606_idmapping_selected.tab"
plot_folder = "../plots"

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
mprotien_go_dict, mgo_protein_dict= parser2.read_in_GO(mouse_go_file,mouse_dict.keys())
#hprotien_go_dict, hgo_protein_dict = parser2.read_in_GO(human_go_file,human_subset)
#reference_protien_go_dict, go_rotein_dict  = parser2.merge_GO(mprotien_go_dict,homolog_dict,hprotien_go_dict)

gomp_args = [ 45.3960489 ,6.75244581 ,1.21480502]
cor_dict = align.error_function_summarize_dict(curves.gom,homolog_dict,mouse_dict,
                                                      align.curve_pearson, times_uniq,
                                                      spline_dict, max,
                                                      weight_dict, gomp_args)

print len(mgo_protein_dict)
print len(mprotien_go_dict)
print cor_dict.keys()[1]
print mgo_protein_dict.items()[1]
go_total_dict,go_thresh_dict = significance.go_analyis(cor_dict, mgo_protein_dict, .7)
print len(go_total_dict)
with open("go.csv",mode='wb') as gofile:
    gowriter = csv.writer(gofile)
    for go in go_total_dict.keys():
        gowriter.writerow([go,go_total_dict[go], go_thresh_dict[go]])