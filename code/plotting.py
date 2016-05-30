import numpy as np
import matplotlib.pyplot as plt

def plot_spline(times, spline_fun, folder, name):
    filename = folder+"/"+name+".png"
    xvals = np.linspace(times[0], times[-1], num=200)
    yvals = [spline_fun(xval) for xval in xvals]
    plt.plot(xvals,yvals)
    plt.xlabel("days")
    plt.ylabel("expression")
    plt.title(name)
    plt.savefig(filename)
    plt.clf()
