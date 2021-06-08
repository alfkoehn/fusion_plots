# coding=utf-8

"""
Plot the plasma zoo diagramm (plasma-types in T-n diagramm).

Simply run this script to produce a png plot:
    $ python plasma_zoo.py
"""

__author__      = 'Alf Köhn-Seemann'
__email__       = 'koehn@igvp.uni-stuttgart.de'
__copyright__   = 'University of Stuttgart'
__license__     = 'MIT'

# import standard modules
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as consts

from matplotlib.colors import LogNorm
from matplotlib import ticker


def calc_debye( n=1e20, T=1, unit='eV' ):
#;{{{
    """
    Calculate the Debye length.

    Parameters
    ----------
    n: float
        plasma density in m^-3
    T: float
        plasma temperature in K (or eV, see parameter 'unit')
    unit: str
        if set to 'eV', plasma temperature is assumed to be in eV

    Returns
    -------
    float
        Debye length in meters.
    """

    if unit == 'eV':
        T *= consts.e/consts.k

    return np.sqrt( consts.epsilon_0 * consts.k * T / (consts.e**2 * n) )
#;}}}


def calc_ND( n=1e20, T=1, unit='eV' ):
#;{{{
    """
    Calculate the plasma parameter (number of particles in Debye sphere).

    Parameters
    ----------
    n: float
        plasma density in m^-3
    T: float
        plasma temperature in K (or eV, see parameter 'unit')
    unit: str
        if set to 'eV', plasma temperature is assumed to be in eV

    Returns
    -------
    float
        Number of particles in Debye sphere.
    """

    lambda_D = calc_debye(n,T,unit=unit)

    return n * 4./3. * np.pi * lambda_D**3

#;}}}


def calc_Trel():
#;{{{
    """
    Calculate the temperature when a plasma becomes relativistic.

    Parameters
    ----------

    Returns
    -------
    float
        Temperature in eV above which the plasma becomes relativitic. 
    """

    return consts.m_e*consts.c**2 / consts.e
#;}}}


def calc_Tdeg( plasma_density ):
#;{{{
    """
    Calculate the plasma temperature, when the plasma becomes degenerated.

    Parameters
    ----------
    plasma_density: float
        plasma density in m^-3

    Returns
    -------
    float
        temperature in eV
    """

    return consts.hbar**2/(2.*consts.m_e) * (3.*np.pi**2*plasma_density)**(2./3.) / consts.e

#;}}}


def calc_Tnonideal( plasma_density ):
#;{{{
    """
    Calculate the plasma temperature, when the plasma becomes non-ideal

    Parameters
    ----------
    plasma_density: float
        plasma density in m^-3

    Returns
    -------
    float
        temperature in eV
    """

    # non-ideal plasmas with strong coupling parameter
    return consts.e**2/(4.*np.pi*consts.epsilon_0) * plasma_density**(1./3.) / consts.e

#;}}}


def build_plasma_zoo():
#;{{{
    """
    Return a dictionary containing the plasma zoo.

    The keys of the dictionary are strings labelling the plasma type.
    For each key, a numpy array with two elements is returned, 
    where the first element corresponds to the plasma density,
    the second to the plasma temperature.

    Parameters
    ----------

    Returns
    -------
    dictionary   
    """

    plasma_zoo = {
            'interstellar\nmedium':     np.array([1e7, .8e0]),      #rather:ne=1e6
            'solar\ncore':              np.array([1e30, 1e3]),      #ok
            'ionosphere':               np.array([1e11, 1e-1]),     #ok
            'flames':                   np.array([1e16, 1e-1]),     #ok
            r'e$^{-}$'+'gas in\nmetals':np.array([1e29, 5e-2]),     #ok
            'solar\nwind':              np.array([1e7,  1e1]),      #ok
#            'interplanetary':          np.array([1e11,1e1]),                  #
            'gas\ndischarge':           np.array([5e16, 1e0]),      #ok
            'lightning':                np.array([1e20, 1e0]),      #ok
            'white\ndwarf':             np.array([1e33, 2e0]),      #ok
            'solar\ncorona':            np.array([1e15, 1e2]),      #ok
            'magnetic\nfusion':         np.array([1e20, 1e4]),      #ok
            'inertial\nfusion':         np.array([1e30, 1e4]),      #300-1000g/cm^3 in burn phase = 1e32 
            'magnetosphere\nof pulsars':np.array([1e10, 1e6]),      #ok
            }

    return plasma_zoo

#;}}}


def str_fmt(x):

    power = int(round(np.log10(x)))

    return r'$N_D=10^{{{0}}}$'.format(power)


def write_plasma_zoo_into_plot( ax, plasma_zoo, 
                                plot__lambda_D=False,
                                silent=True 
                              ):
#;{{{
    """
    Write plasma zoo into plot.

    Parameters
    ----------
    ax: Axes object
        Axes object into which the plasma zoo will be written
    plasma_zoo: dict
        dict object which contains the plasma zoo, expected to consist
        of a key (well, obviously, otherwise no dict...) and a two-element
        numpy array: [plasma density in m^-3, plasma temperature in eV]
    plot__lambda_D: bool
        set if also the Debye-length is included in the Axes object,
        a few formatting things will be different then
    silent: bool
        if False, some useful (?) output will be printed to console

    Returns
    -------
    """

    if not silent:
        print( 'write_plasma_zoo_into_plot' )

    if plot__lambda_D:
        plasma_zoo_col = 'dimgrey'
    else:
        plasma_zoo_col = 'black'

    for key in plasma_zoo.keys():
        if not silent:
            print( '    {0}: {1:8.2e} m^-3, {2:8.2e} eV'.format(
                    key.replace('\n', ' '), plasma_zoo[key][0], plasma_zoo[key][1]) )

        ax.text( plasma_zoo[key][0], plasma_zoo[key][1],
                 key,
                 color=plasma_zoo_col,
                 horizontalalignment='center', verticalalignment='center'
               )

#;}}}


def make_lambda_D_contours( fig, ax, 
                            T_vals=[], n_vals=[],
                            silent=True,
                          ):
#;{{{
    """
    Plot filled contours of Debye length into plot. 

    Parameters
    ----------
    fig: Figure object
        Figure object belonging to 'ax' (see parameter below)
    ax: Axes object
        Axes object into which the plasma zoo will be written
    T_vals: numpy array of floats
        plasma temperature in eV, corresponding to y-axis
    n_vals: numpy array of floats
        plasma density in m^-3, corresponding to x-axis
    silent: bool
        if False, some useful (?) output will be printed to console

    Returns
    -------
    """

    fct_name = 'make_lambda_D_contours'

    # check if temperature and density arrays were provided
    # if not, create them
    if len(T_vals) == 0:
        # plasma temperature in eV
        T_vals = np.logspace( np.log10(1e-2), np.log10(1e7),  num=1000 )
    if len(n_vals) == 0:
        # plasma density in m^-3
        n_vals = np.logspace( np.log10(1e5),  np.log10(1e35), num=2000 )

    # spatial coordinates (2D) for contour plot
    nn, TT = np.meshgrid( n_vals, T_vals )

    # caclulate the Debye length
    lambda_D = np.empty( (T_vals.shape[0], n_vals.shape[0] ) )
    for ii in np.arange(n_vals.shape[0]):
        for jj in np.arange(T_vals.shape[0]):
#            print( 'ii={0:d}, jj={1:d}, n={2:13.6e}, T={3:13.6e}, lambda_D={4:13.6e}'.
#                    format( ii, jj, n_vals[ii], T_vals[jj], calc_debye(n=n_vals[ii],T=T_vals[jj]) )
#                 )
            lambda_D[jj,ii] = calc_debye( n=n_vals[ii], T=T_vals[jj] )

    # identify non-ideal plasma
    # relativistic plasmas
    T_rel = calc_Trel()
    # degenerated plasmas
    TT_deg = calc_Tdeg( nn )
    # non-ideal plasmas with strong coupling parameter
    T_nonideal = calc_Tnonideal( nn )
    # get indices of non-ideal plasmas in spatial coordinates 
    TT_rel_ids      = (TT >= T_rel)
    TT_deg_ids      = (TT <= TT_deg)
    TT_nonideal_ids = (TT <= T_nonideal)

    # set lambda_D at non-ideal plasma to NaN in order to not plot it 
    lambda_D[TT_rel_ids]      = np.nan
    lambda_D[TT_deg_ids]      = np.nan
    lambda_D[TT_nonideal_ids] = np.nan

    # contour levels are logarithmic due to large range
    lD_contLevels = np.logspace( np.log10(1e-12), 
                                 np.log10(1e4), 
                                 9 )

    if not silent:
        print( '{0}: lambda_D contour levels:'.format(fct_name) )
        print( lD_contLevels )

    cont_lD = ax.contourf( nn, TT, lambda_D,
                           levels=lD_contLevels,
                           norm=LogNorm()
                         )
    locator = ticker.LogLocator(base=10)

    # add colorbar
    cbar = fig.colorbar( cont_lD, fraction=0.046, pad=0.04, ticks=locator )
    cbar.ax.tick_params( direction='in' )
    cbar.set_label( 'Debye length in m' )

#;}}}


def make_N_D_contours( fig, ax, 
                       T_vals=[], n_vals=[],
                       silent=True,
                     ):
#;{{{
    """
    Plot contour levels of plasma parameter into plot. 

    Parameters
    ----------
    fig: Figure object
        Figure object belonging to 'ax' (see parameter below)
    ax: Axes object
        Axes object into which the plasma zoo will be written
    T_vals: numpy array of floats
        plasma temperature in eV, corresponding to y-axis
    n_vals: numpy array of floats
        plasma density in m^-3, corresponding to x-axis
    silent: bool
        if False, some useful (?) output will be printed to console

    Returns
    -------
    """

    fct_name = 'make_N_D_contours'

    # check if temperature and density arrays were provided
    # if not, create them
    if len(T_vals) == 0:
        # plasma temperature in eV
        T_vals = np.logspace( np.log10(1e-2), np.log10(1e7),  num=1000 )
    if len(n_vals) == 0:
        # plasma density in m^-3
        n_vals = np.logspace( np.log10(1e5),  np.log10(1e35), num=2000 )

    # spatial coordinates (2D) for contour plot
    nn, TT = np.meshgrid( n_vals, T_vals )

    # calculate plasma parameter
    N_D = np.empty( (T_vals.shape[0], n_vals.shape[0] ) )
    for ii in np.arange(n_vals.shape[0]):
        for jj in np.arange(T_vals.shape[0]):
            N_D[jj,ii] = calc_ND( n=n_vals[ii], T=T_vals[jj] )

    # identify non-ideal plasma
    # relativistic plasmas
    T_rel = calc_Trel()
    # degenerated plasmas
    TT_deg = calc_Tdeg( nn )
    # non-ideal plasmas with strong coupling parameter
    T_nonideal = calc_Tnonideal( nn )
    # get indices of non-ideal plasmas in spatial coordinates 
    TT_rel_ids      = (TT >= T_rel)
    TT_deg_ids      = (TT <= TT_deg)
    TT_nonideal_ids = (TT <= T_nonideal)

    # set N_D at non-ideal plasma to NaN in order to not plot it 
    N_D[TT_rel_ids]      = np.nan
    N_D[TT_deg_ids]      = np.nan
    N_D[TT_nonideal_ids] = np.nan

    # contour levels are logarithmic due to large range covered
    ND_contLevels = np.logspace( np.log10(1e0), 
                                 np.log10(1e15), 
                                 6 )

    if not silent:
        print( '{0}: N_D contour levels:'.format(fct_name) )
        print( ND_contLevels )

    # manually set position for labels of contour levels
    ND_contLabelsPos = [ (1e26,3.5e1), 
                         (1e31,2e5),
                         (1e25,2e5),
                         (1e19,2e5),
                         (1e13,2e5),
                         (1e7, 2e5)
                       ]
    cont_ND = ax.contour( nn, TT, N_D,
                          levels=ND_contLevels,
                          colors='darkgrey', linestyles='dashed',
                        )

    # NOTE: EVIL HACK to manually write contour label
    #       reason was that clabels was not working properly
    #       probably due to setting some areas to NaN
    for ii in np.arange(len(ND_contLabelsPos)):
        ax.text( ND_contLabelsPos[ii][0], ND_contLabelsPos[ii][1], 
                 str_fmt(ND_contLevels[ii]),
                 rotation=40, 
                 fontsize=10, color='darkgrey'
               )
        if not silent:
            print( '{0}: {1}, contour-level = {2}, formatted string-label = {3}'.format(
                    fct_name, ii, ND_contLevels[ii], str_fmt(ND_contLevels[ii])) )

#;}}}


def write_plasma_limits_into_plot( ax, 
                                   plot__lambda_D=False, xkcd_style=True, 
                                   T_vals=[], n_vals=[],
                                   silent=True
                                 ):
#;{{{
    """
    Mark (and label) limit of ideal plasmas in plot. 

    Parameters
    ----------
    ax: Axes object
        Axes object into which the plasma zoo will be written
    plot__lambda_D: bool
        set if also the Debye-length is included in the Axes object,
        a few formatting things will be different then
    xkcd_style: bool
        set True, if xkcd plot style is used,
        a few formatting things will be different then
    T_vals: numpy array of floats
        plasma temperature in eV, corresponding to y-axis
    n_vals: numpy array of floats
        plasma density in m^-3, corresponding to x-axis
    silent: bool
        if False, some useful (?) output will be printed to console

    Returns
    -------
    """

    fct_name = 'write_plasma_limits_into_plot'

    # check if temperature and density arrays were provided
    # if not, create them
    if len(T_vals) == 0:
        # plasma temperature in eV
        T_vals = np.logspace( np.log10(1e-2), np.log10(1e7),  num=1000 )
    if len(n_vals) == 0:
        # plasma density in m^-3
        n_vals = np.logspace( np.log10(1e5),  np.log10(1e35), num=2000 )

    # spatial coordinates (2D) for contour plot
    nn, TT = np.meshgrid( n_vals, T_vals )
 
    # label boundary for relativistic plasmas
    ax.hlines( y=calc_Trel(), xmin=np.nanmin(nn), xmax=np.nanmax(nn),
               linestyles='solid', linewidth=3, colors='grey' )
    ax.text( 1e20, 9e5, 'relativistic plasmas', color='grey' )

    # label boundary for degenerated plasmas
    ax.plot( n_vals, calc_Tdeg(n_vals),
             linestyle='solid', linewidth=3, color='grey' )
    label_deg_n = 5e30
    label_deg_T = 2e4
    # failed attemp to make rotation fit to T_deg-function
    label_deg_n_id = np.where( np.abs(n_vals-label_deg_n) == np.abs(n_vals-label_deg_n).min() )
    label_deg_n_id = label_deg_n_id[0][0]
    T_deg_vals = calc_Tdeg(n_vals)
    ## angle in data coordinates
    label_deg_angle_data = np.rad2deg( np.arctan2( T_deg_vals[label_deg_n_id] - T_deg_vals[(label_deg_n_id-1)],
                                                   n_vals[label_deg_n_id]     - n_vals[(label_deg_n_id-1)]) )
    ## angle in screen coordinates
    label_deg_angle_screen = ax.transData.transform_angles( np.array((label_deg_angle_data,)),
                                                            np.array([n_vals[(label_deg_n_id-1)],
                                                                      T_deg_vals[(label_deg_n_id-1)]]).reshape((1,2)))[0]
#    ax.annotate( 'Text', 
#                 xy=(n_vals[label_deg_n_id-1], T_deg_vals[label_deg_n_id-1]),
#                 rotation_mode='anchor', rotation=label_deg_angle_screen
#               )
    if not silent:
        print( fct_name )
        print( label_deg_n, label_deg_n_id, n_vals[label_deg_n_id], T_deg_vals[label_deg_n_id] )
        print( n_vals[label_deg_n_id-1], T_deg_vals[label_deg_n_id-1] )
        print( label_deg_angle_data, label_deg_angle_screen )
    if plot__lambda_D:
        label_deg_angle      = 62.
        label_nonideal_angle = 42.
    else:
        label_deg_angle      = 59.
        label_nonideal_angle = 39.
    if xkcd_style:
        label_deg_n = 3e29
        label_nonideal_T = 1.5e0
    else:
        label_nonideal_T = 5e-1
    label_nonideal_n = 3e21
    ax.text( label_deg_n, label_deg_T, 
             'degenerated plasmas', 
             rotation=label_deg_angle,
             color='grey' )

    # label boundary to non-ideal plasmas with strong coupling
    ax.plot( n_vals, calc_Tnonideal( n_vals ), 
             linestyle='solid', linewidth=3, color='grey' )
    ax.text( label_nonideal_n, label_nonideal_T, 
             'non-ideal plasmas', 
             rotation=label_nonideal_angle,
             color='grey' )

#;}}}


def main():
#;{{{

    print( '\n' )
    print( 'Let me know if you have questions, requests or found some bugs.')
    print( '                                    -- Alf Köhn-Seemann, April 2020\n' )

    plot__lambda_D = True
    plot__N_D      = True
    plot__limits   = True
    label_plasmas  = True

    # plasma temperature in eV, plasma density in m^-3
    T_vals = np.logspace( np.log10(1e-2), np.log10(1e7),  num=1000 )
    n_vals = np.logspace( np.log10(1e5),  np.log10(1e35), num=2000 )

    fname_plot = ''
    xkcd_style = True

    # plot configuration
    # optionally acitvate xkcd-style plot
    if xkcd_style:
        plt.xkcd()

    fig1 = plt.figure( figsize=(8,6) )
    ax1  = fig1.add_subplot( 1,1,1 )

    if plot__lambda_D:
        make_lambda_D_contours( fig1, ax1, 
                                T_vals=T_vals, n_vals=n_vals,
                                silent=True,
                              )

    if plot__N_D:
        make_N_D_contours( fig1, ax1, 
                           T_vals=T_vals, n_vals=n_vals,
                           silent=True,
                         )

    if label_plasmas:
        # get the plasma zoo
        plasma_zoo = build_plasma_zoo()

        # for xkcd-style, a small correction is necessary
        # otherwise, the following label would overlap with another
        # (due to different font size)
        if xkcd_style:
            plasma_zoo['lightning'][0] = 5e21
    
        write_plasma_zoo_into_plot( ax1, plasma_zoo, plot__lambda_D )

    if plot__limits:
        write_plasma_limits_into_plot( ax1, 
                                       plot__lambda_D=plot__lambda_D, xkcd_style=xkcd_style, 
                                       T_vals=T_vals, n_vals=n_vals,
                                     )

    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax1.set_xlim( np.nanmin(n_vals), np.nanmax(n_vals) )
    ax1.set_ylim( np.nanmin(T_vals), np.nanmax(T_vals) )

    ax1.set_xticks([1e5,1e10,1e15,1e20,1e25,1e30,1e35])
    ax1.set_yticks([1e-2,1e0,1e2,1e4,1e6])

    ax1.set_xlabel( r'plasma density in m$^{-3}$' )
    ax1.set_ylabel( r'temperature in eV' )

    # force ticks to point inwards
    ax1.tick_params( axis='both', which='both', direction='in',
                     top=True, right=True
                   )

    ax1.minorticks_off()

    # write credits into plot, to ensure that people know they can use the plot
    # (somebody once told me, every plot appearing somewhere in the internet
    #  should contain information on how to use it, otherwise it is useless)
    # you probably want to remove it when you make you own plot
    # attribution would still be gratefully acknowledged :)
    # also note that the licence refers only to that specific plot
    # the licence for the code is mentioned above and in the LICENCE file
    credit_str = u'Alf Köhn-Seemann, CC BY-SA 4.0'
    fig1.text( .65, .89, credit_str, fontsize=7 )

    if len( fname_plot ):
        plt.savefig( fname_plot, bbox_inches='tight', dpi=600 )
        print( '    plot written into {0}'.format( fname_plot ) )
    else:
        plt.show()

#;}}}


if __name__ == '__main__':
    main()

