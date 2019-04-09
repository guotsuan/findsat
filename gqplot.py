#import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoLocator, AutoMinorLocator, MultipleLocator, FormatStrFormatter


def set_plot(ax, label_size=None):
    """ gqplot module
    set the plot parameters for paper plot """

    #ax.xticks(size=18)
    #ax.yticks(size=18)
    if label_size is None:
        ax.tick_params(axis='both', labelsize=16)
    else:
        ax.tick_params(axis='both', labelsize=label_size)
    ax.tick_params(axis='both', pad=10.0)

    ax.spines['left'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)
    ax.tick_params(axis='both', width=2.0)
    #ax.tick_params(axis='both', which='major', length=10)
    #ax.tick_params(axis='both', which='minor', length=5)
    ax.tick_params(axis='both', which='major', length=6)
    ax.tick_params(axis='both', which='minor', length=3)
    ax.tick_params(which='minor', width=1.5)
    ax.xaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ylim_low, ylim_up = ax.get_ylim()
    if np.abs(ylim_low - ylim_up) < 1.0:
        ax.yaxis.set_major_locator(AutoLocator())
        myfmt = FormatStrFormatter('%3.2f')
        ax.yaxis.set_major_formatter(myfmt)
    else:
        ax.yaxis.set_major_locator(AutoLocator())
        #ax.yaxis.set_major_locator(MultipleLocator(1))



def logerr(y, err):

    """ logerr
    y: y coordinate data
    err: err of y data """

    low_err = np.log10(y) - np.log10(np.clip(y - err, 1e-5, 1e10))
    up_err = np.log10(y + err) - np.log10(y)
    return [low_err, up_err]


    #allerr=np.zeros((2,mock_3dx.shape[0]))+1.0
    #allerr[1,:]=np.log10(mock_3dy+mock_3derr)-np.log10(mock_3dy)
    #allerr[0,:]=np.log10(mock_3dy)-np.log10(mock_3dy-mock_3derr)
def gama_z00_phi(x):
    m_star = 20.73
    phi_star = 0.90
    alpha = -1.22
    result = 0.0
    result = 0.4 * np.log(10) * phi_star * ((10. ** (m_star - x))) ** (1 + alpha) * \
        np.exp(-10 ** (0.4 * (m_star - x)))
    return result
