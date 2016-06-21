import align
from multiprocessing import Lock,pool
import random
from scipy import optimize
import math
def parallel_random_comp(funct,init,lock,filename,seed):
    options = {"method":"Nelder-Mead"}
    seed_funct = lambda x: funct(seed, x)
    result = optimize.basinhopping(seed_funct, init, minimizer_kwargs=options)
    lock.acquire(True)
    with open(filename, 'a') as write_file:
        write_file.write(str(seed_funct(result.x))+"\n")
    lock.release()

def FDR(funct,init,filename,times):
    write_lock = Lock()
    random_seeds = [random.random() for i in xrange(times)]
    print random_seeds[1:5]
    par_comp = lambda seed: parallel_random_comp(funct,init,write_lock, filename,seed)
    tpool = pool.ThreadPool()
    tpool.map(par_comp,random_seeds)





def go_analyis(cor_dict, go_protein_dict, thresh):
    go_total_dict = {}
    go_thresh_dict = {}
    for go_cat in go_protein_dict.keys():
        go_total_dict[go_cat] = 0
        go_thresh_dict[go_cat] = 0
        for protein in go_protein_dict[go_cat]:
            if protein in cor_dict:
                go_total_dict[go_cat] = go_total_dict[go_cat]+1
                if cor_dict[protein] >= thresh:
                    go_thresh_dict[go_cat] = go_thresh_dict[go_cat]+1
    return go_total_dict,go_thresh_dict