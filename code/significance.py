import align

def FDR(human_proteins, htimes, error_fun, opt_fun, init_param):
    vals = []
    for i in xrange(1000):
        spline_dict = align.make_spline_dict_rand(human_proteins, htimes)
        goal_fun = lambda x: error_fun(spline_dict,x)
        result = opt_fun(goal_fun, init_param)
        vals.append(goal_fun(result.x))
    return vals

def find_best(error_function, init_param):
    return None