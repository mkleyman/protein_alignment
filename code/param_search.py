import numpy as np

def record_twoparam(eval_fun, range1, range2):
    record_arr = np.zeros(shape=(len(range1)+1,len(range2)+1))
    print record_arr.shape
    x = 0
    for i in range1:
        y = 0
        for j in range2:
            y += 1
            record_arr[x,y] =  eval_fun([i,j])
        x+= 1
    return record_arr