#!/usr/bin/env python
# coding=utf-8

"""
Plot the triple product as a function of time.
"""

__author__      = 'Alf Köhn-Seemann'
__email__       = 'koehn@igvp.uni-stuttgart.de'
__copyright__   = 'Uni Stuttgart'
__license__     = 'MIT'

# import standard modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy.polynomial.polynomial as poly

# write credits into plot, to ensure that people know then can use the plot
# I once learned that every plot appearing somewhere on the internet should
# contain infos how to use, otherwise it is useless.
# Note that the license refers only to the plot, the license for the code
# is mentioned above and in the LICENSE file
credit_str  = f'{__author__}, CC BY-SA 4.0'


def get_dataset( dataset='Webster' ):
#{{{
    """
    Return nTtau vs time dataset.

    Ikeda paper: Ikeda, NF 50 (2010) paper, Fig. 1
        https://doi.org/10.1088/0029-5515/50/1/014002
    Webster-paper: Webster, Phys. Educ. (2003), Fig. 3
        https://doi.org/10.1088/0031-9120/38/2/305
    Dataset extracted using the WebPlotDigitizer.
    """

    possible_datasets   = ['Ikeda', 'Webster']

    if dataset == 'Ikeda':
        data = [
            [ 1968, 1.20724640e-3, 'T3' ],
            [ 1971, 5.64724637e-3, 'ST' ],
            [ 1975, 2.22299648e-2, 'TFR'], 
            [ 1978, 6.19673987e-2, 'PLT'], 
            [ 1978, 2.43930051e-1, 'Alcator A'], 
            [ 1981, 1.11644346e-1, 'PDX'], 
            [ 1983, 1.21270482e+0, 'Alcator C'], 
            [ 1984, 8.58771357e-1, 'JET'], 
            [ 1984, 4.62359294e-1, 'DIII'], 
            [ 1986, 1.20045791e+0, 'JET'], 
            [ 1986, 1.82004293e+0, 'TFTR'], 
            [ 1989, 9.14070287e+0, 'JT60U'], 
            [ 1991, 1.29079461e+1, 'JT60U'], 
            [ 1992, 4.40803362e+1, 'JT60U'], 
            [ 1993, 1.06599484e+2, 'JT60U'], 
            [ 1994, 2.14742679e+2, 'JT60U'], 
            [ 1996, 1.46018594e+2, 'JT60U'], 
            [ 1996, 8.02281043e+1, 'TFTR'], 
            [ 1996, 5.97709001e+1, 'DIII-D'], 
            [ 1998, 1.16796162e+2, 'JET']
            ]
    elif dataset == 'Webster':
        data = [
            [ 1968, 0.0011831415917701873, 'T3'], 
            [ 1971, 0.005544496169709835,  'ST'], 
            [ 1975, 0.022003379868766233, 'TFR'], 
            [ 1978, 0.06186522810656025,  'PLT'], 
            [ 1978, 0.24520180197227637,  'Alcator A'], 
            [ 1981, 0.11280221433667462,  'PDX'], 
            [ 1983, 1.2020523507668979,   'Alcator C'], 
            [ 1984, 0.45755718141889035,  'DIII'], 
            [ 1984, 0.8427948172631539,   'JET'], 
            [ 1986, 1.1899113127482666,   'JET'], 
            [ 1986, 1.7946141450897937,   'TFTR'], 
            [ 1989, 9.191399444917563,    'JT-60U'], 
            [ 1991, 12.977200281296101,   'JT-60U'], 
            [ 1992, 45.03113151570462,    'JT-60U'], 
            [ 1993, 107.11518319232955,   'JT-60U'], 
            [ 1994, 218.1080833047305,    'JT-60U'], 
            [ 1996, 62.220159508278485,   'DIII-D'], 
            [ 1996, 80.32552821013279,    'TFTR'], 
            [ 1996, 149.6097332253447,    'JT-60U'], 
            [ 1998, 117.25347357319447,   'JET']
            ]
    else:
        print( 'ERROR: dataset does not exist' )
        print( '       possible datasets are ', possible_datasets )

    return data
#}}}


def extract_data( dataset, data2extract ):
#{{{
    """
    Extracts the desired data from one of the datasets. Possible value are 
    'year', 'nTtau', 'name'.
    """

    possible_data2extract   = ['year', 'nTtau', 'name']

    if data2extract == 'year':
        chNr    = 0
    elif data2extract == 'nTtau':
        chNr    = 1
    elif data2extract == 'name':
        chNr    = 2
    else:
        print( 'ERROR: only the following data can be extracted: ')
        print( '       ', possible_data2extract )
        return -1

    # more efficient to first copy every value into a list and then transfort
    # into numpy array, instead of using numpy functions like append or
    # vstack or hstack (because everytime one of those numpy functions is
    # called, the data is copied)
    data2return = []

    for ii in range(len(dataset)):
        data2return.append( dataset[ii][chNr] )

    return( np.array(data2return) )

#}}}


def make_plot( fname_plot='' ):
#{{{
    """
    Output a plot, either to X-window (default) or into file.

    Parameters
    ----------
    fname_plot: str
        if empty, plot will be output into X-window.

    Returns
    -------
    """

    if len(fname_plot) > 0:
        plt.savefig( fname_plot, dpi=600, bbox_inches='tight' )
        print( 'written plot into file {0}'.format(fname_plot) )
    else:
        plt.show()
#}}}


def plot_nTtau_time( dataset='Webster', add_ITER=True, make_fit=True, 
                     fname_plot='' ):
#{{{
    """
    """

    # get the data and scale it where necessary
    data        = get_dataset( dataset=dataset )
    name_vals   = extract_data( data, 'name' )
    year_vals   = extract_data( data, 'year' )
    nTtau_vals  = extract_data( data, 'nTtau' )

    # prepare plot
    fig = plt.figure( figsize=(8,6) )       # (width, height) in inches
    ax1 = fig.add_subplot( 1,1,1 )          # rows, cols, plotnumber

    # optionally, perform and plot a linear fit to log(nTtau) dataset
    # do this without the ITER datapoint
    if make_fit:
        # perform a linear fit to the log(nTtau) data
        x_new = np.linspace( np.min(year_vals), np.max(year_vals), 100 )
        coefs = poly.polyfit( year_vals, np.log10(nTtau_vals), 1 )
        fit   = poly.polyval( x_new, coefs )
        # calculate the doubling time
        t_double = np.log10(2)/coefs[1]
        print( 'fit performed to y  =  a0 + a1*x, a0={}, a1={}'.format( coefs[0], coefs[1] ) )
        print( '    where  y  --> log10(nTtau)' )
        print( '           a0 --> log10(nTtau)_0' )
        print( '           a1 --> log10(a)' )
        print( '           x  --> t' )
        print( '    and  (nTtau) = (nTtau)_0 * a^t ' )
        print( '         --> exponential growth' )
        print( '    doubling time: 2*y_0 = y0*a^t ==> 2=a^t' )
        print( '                                  ==> t=log(2)/log(a)={}'.format( t_double ))
        print( '                   --> {} years'.format( np.round( t_double, decimals=1 )))

        # plot fit to data
        ax1.plot( x_new, 10**(fit) )
        # annotate fit
        ax1.annotate( '{0:3.1f} years doubling rate (fit to data)'.format(np.round(t_double, decimals=1)), 
                      xy=( x_new[20], 10**(fit[20]) ), xytext=( 1976, .6*nTtau_vals[1] ),
                      arrowprops=dict( arrowstyle="->", shrinkA=0, color='0.2' ), 
                      ha='left', va='top'
                    )

    # optionally, add ITER point
    if add_ITER:
        nTtau_vals  = np.append(nTtau_vals, 4e21*1e-19)
        year_vals   = np.append(year_vals, 2035)
        name_vals   = np.append(name_vals, 'ITER')

    # plot nTtau as a function time dataset
    ax1.plot( year_vals, nTtau_vals, 'o', markersize=10 )

    # plot format stuff
    ax1.set_xlabel( "Year" , fontsize=16 )
    ax1.set_ylabel( r"$nT_i\tau_E$ in ($10^{19}\,$keVm$^{-3}$s)", fontsize=16 )
    ax1.set_ylim( 1e-3, 1e3 )
    ax1.set_yscale('log')
    ax1.xaxis.set_major_formatter( FormatStrFormatter('%4d') )
    ax1.tick_params( axis='both', which='both', direction='in', 
                     labelsize=14, 
                     top='on', right='on' )

    # add annotations to datapoints in a very unefficient (annoying) way
    # changing the DPI of the image requires to adjust the pixel-offset values
    # this is realized by the factor dpi_scale, 1 corresponds to dpi=100
    if len(fname_plot) == 0:
        dpi_scale = 1
    else:
        dpi_scale = 5
    # loop through dataset
    for ii in range( name_vals.shape[0] ):
        # define list for cases which cannot be annotated automatically
        extra_labels = [ 'Alcator A', 'Alcator C', 'JT60U', 'JT-60U', 'JET', 'ITER' ]
        # automatic annotation
        if name_vals[ii] not in extra_labels:
            if add_ITER:
                ax1.annotate( name_vals[ii], 
                              xy=( year_vals[ii], nTtau_vals[ii]), 
                              xytext=(10*dpi_scale,-5*dpi_scale), textcoords='offset pixels',
                            )
            else:
                ax1.annotate( name_vals[ii], 
                              xy=( year_vals[ii], nTtau_vals[ii]), 
                              xytext=(10*dpi_scale,-5*dpi_scale), textcoords='offset pixels',
                            )
        # handles the above defined extra cases
        elif name_vals[ii] == 'Alcator A' or name_vals[ii] == 'Alcator C': 
            if add_ITER:
                ax1.annotate( name_vals[ii], 
                              xy=( year_vals[ii], nTtau_vals[ii]), 
                              xytext=(-85*dpi_scale,-6*dpi_scale), textcoords='offset pixels',
                            )
            else:
                ax1.annotate( name_vals[ii], 
                              xy=( year_vals[ii], nTtau_vals[ii]), 
                              xytext=(-90*dpi_scale,-6*dpi_scale), textcoords='offset pixels',
                            )
        elif (name_vals[ii] == 'JT60U') or (name_vals[ii] == 'JT-60U'):
            # write left of symbol
            if (year_vals[ii] < 1990) or ( year_vals[ii] > 1991 and year_vals[ii] < 1995):
                if add_ITER:
                    ax1.annotate( name_vals[ii], 
                                  xy=( year_vals[ii], nTtau_vals[ii]), 
                                  xytext=(-65*dpi_scale,-6*dpi_scale), textcoords='offset pixels',
                                )
                else:
                    ax1.annotate( name_vals[ii], 
                                  xy=( year_vals[ii], nTtau_vals[ii]), 
                                  xytext=(-64*dpi_scale,-6*dpi_scale), textcoords='offset pixels',
                                )
            # label above symbol
            elif year_vals[ii] > 1995:
                ax1.annotate( name_vals[ii], 
                              xy=( year_vals[ii], nTtau_vals[ii]), 
                              xytext=(-20*dpi_scale,7*dpi_scale), textcoords='offset pixels',
                            )
            # default case (label right of symbol)
            else:
                if add_ITER:
                    ax1.annotate( name_vals[ii], 
                                  xy=( year_vals[ii], nTtau_vals[ii]), 
                                  xytext=(12*dpi_scale,-5*dpi_scale), textcoords='offset pixels',
                                )
                else:
                    ax1.annotate( name_vals[ii], 
                                  xy=( year_vals[ii], nTtau_vals[ii]), 
                                  xytext=(10*dpi_scale,-5*dpi_scale), textcoords='offset pixels',
                                )
        elif name_vals[ii] == 'JET':
            if year_vals[ii] > 1995:
                ax1.annotate( name_vals[ii], 
                              xy=( year_vals[ii], nTtau_vals[ii]), 
                              xytext=(-6*dpi_scale,9*dpi_scale), textcoords='offset pixels',
                            )
            else:
                if add_ITER:
                    ax1.annotate( name_vals[ii], 
                                  xy=( year_vals[ii], nTtau_vals[ii]), 
                                  xytext=(12*dpi_scale,-5*dpi_scale), textcoords='offset pixels',
                                )
                else:
                    ax1.annotate( name_vals[ii], 
                                  xy=( year_vals[ii], nTtau_vals[ii]), 
                                  xytext=(10*dpi_scale,-5*dpi_scale), textcoords='offset pixels',
                                )
        elif name_vals[ii] == 'ITER':
            ax1.annotate( name_vals[ii], 
                          xy=( year_vals[ii], nTtau_vals[ii]), 
                          xytext=(-48*dpi_scale,-6*dpi_scale), textcoords='offset pixels',
                        )

    fig.text( .703, .885, credit_str, fontsize=7 )

    make_plot( fname_plot )

#}}}


def main():

    plot_nTtau_time( add_ITER=True, 
                     fname_plot='nTtau_vs_time.png'
                   )


if __name__ == '__main__':
    main()

