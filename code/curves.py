import math
import numpy as np

def linear(x,coef):
    return float(coef[1])*float(x)+float(coef[0])


def logit(x,coef):
    return coef[0]*math.log1p(coef[1]/(coef[2]-coef[3]*x))+coef[4]

def exp(x,coef):
    return coef[0]*math.pow(math.e,coef[2]*x)+coef[1]

def exp4(x,coef):
    return coef[0] * math.pow(math.e, coef[2] * x+ coef[3]) + coef[1]
def exp2(x,coef):
    return math.pow(math.e,coef[1]*x)+coef[0]

def quadratic(x,coef):
    return float(coef[1])*float(x)+float(coef[0])+float(coef[2])*math.pow(x,2.0)

def cubic(x,coef):
    return float(coef[1]) * float(x) + float(coef[0]) +\
           float(coef[2]) * math.pow(x, 2.0) + float(coef[3]) * math.pow(x, 3.0)

def cubic2(x,coef):
    return float(coef[1]) * math.pow(x, 3.0) + float(coef[0])

def square_root(x,coef):
    return float(coef[1])*math.sqrt(float(x)+float(coef[2]))+float(coef[0])

def logarithm(x,coef):
    return float(coef[1])*math.log1p(float(coef[2])*x+float(coef[3]))+float(coef[0])