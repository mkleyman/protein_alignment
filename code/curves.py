import math
def linear(x,coef):
    return float(coef[1])*float(x)+float(coef[0])


def logit(x,coef):
    return coef[0]*math.log1p(coef[1]/(coef[2]-coef[3]*x))+coef[4]

def exp(x,coef):
    return coef[0]*math.pow(math.e,coef[2]*x)+coef[1]

def exp2(x,coef):
    return math.pow(math.e,coef[1]*x)+coef[0]
