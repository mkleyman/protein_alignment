class Pars:
    def __init__(self,funct,
               homolog_dict, mouse_dict,
               dif_function, times,
               spline_dict, summary_fun, weight_dict):
        self.funct = funct
        self.homolog_dict = homolog_dict
        self.reference_dict = mouse_dict
        self.dif_function = dif_function
        self.times = times
        self.spline_dict = spline_dict
        self.summary_fun = summary_fun
        self.weight_dict = weight_dict

