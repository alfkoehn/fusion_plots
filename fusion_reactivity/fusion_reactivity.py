# coding=utf-8

"""
Plot the fusion reactivity <sigma*v> as a function of T_ion for various fusion reactios.

Simply run this script to produce a png plot:
    $ python fusion_reactivity.py
"""

__author__      = 'Alf Köhn-Seemann'
__email__       = 'koehn@igvp.uni-stuttgart.de'
__copyright__   = 'University of Stuttgart'
__license__     = 'MIT'

# import standard modules
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as consts
import scipy.interpolate as interp

# change some default properties of matplotlib
plt.rcParams.update( {'font.size':12} )
# force ticks to point inwards
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top']       = True
plt.rcParams['ytick.right']     = True


def make_plot( fname_plot='' ):
#;{{{
    '''
    Output a plot, either to X-window (default) or into file. 

    Parameters
    ----------
    fname_plot: str
        if empty, plot will be output into X-window

    Returns
    -------
    '''


    if len(fname_plot) > 0:
        plt.savefig( fname_plot, dpi=600, bbox_inches='tight' )
        print( 'written plot into file {0}'.format(fname_plot) )
    else:
        plt.show()

#;}}}


def reaction_int2str( reaction_int, silent=True ):
#;{{{
    """
    Translate reaction number to reaction string.

    The idea of this functions is to ensure that in all function within
    this file the same reaction is meant when referring to reaction 3
    (or any other number).

    The integer provided to this function is tranlated into a string
    specifying the reaction:
        1: D + T     --> n + 4He     T(d,n)4He
        2: D + D     --> p + T       D(d,p)T
        3: D + D     --> n + 3He     D(d,n)3He
        4: 3He + D   --> p + 4He     3He(d,p)4He
        5: T + T     --> 4He + n + n
        6: T + 3He   --> D + 4He         41 %
                         p + 4He + n     55 %
                         p + 4He + n     4 %
        7: 3He + 3He --> p + p + 4He 
        8: p + 11B   --> 4He + 4He + 4He
        9: p + p     --> 2H + nu + e^+

    Parameters
    ----------
    reaction_int: int
        possible values see above

    Returns
    -------
    reaction_str: str
        possible values see the dictionary values
    """

    func_name = 'reaction_int'

    reaction_dict = { 1 : 'DT'  , 
                      2 : 'DD_a',
                      3 : 'DD_b',
                      4 : '3HeD',
                      5 : 'TT',
                      6 : 'T3He',
                      7 : '3He3He',
                      8 : 'p11B',
                      9 : 'pp',
                     }

    # check if reaction_int is in keys of dict, return NaN if not
    if reaction_int in reaction_dict.keys():
        reaction_str = reaction_dict[ reaction_int ]
    else:
        reaction_str = 'NaN'
        print( '{0}: ERROR, no reaction for found for provided key'.format( func_name ) )
        print( '{0}  key was {1}'.format( ' '*len(func_name), reaction_int ) )

    if not silent:
        print ('{0}: reaction_int = {1:d}, reaction_str = {2}'.format(
                func_name, reaction_int, reaction_str))

    return reaction_str

#;}}}


def get_fusion_reactivity_Hively( T_ion, reaction=1, extrapolate=True, silent=True ):
#;{{{
    '''
    Return the fusion reactivity for various fusion reactions. 

    A fit to tabulated experimental values is used to get a functional
    dependence of the fusion reactivity on the ion temperature. Fit
    parameters from the following paper are used:
        L.M. Hively, Nuclear Fusion, Vol. 17, No. 4 (1977)
        https://doi.org/10.1088/0029-5515/17/4/019
    Equation (5), referred to as S_5 in the paper, is used here, 
    as it resulted in the best fit.
    Curve fitting range in the paper: T_ion = 1...80 keV
        ==> T_ion values outside of this range should not be used

    Parameters
    ----------
    T_ion: float
        ion temperature in keV, valid range 1-80 keV
    reaction: int
        defines reaction considered, default value is 1
        possible values are:
        1: T + D    --> n + 4He     T(d,n)4He
        2: D + D    --> p + T       D(d,p)T
        3: D + D    --> n + 3He     D(d,n)3He
        4: 3He + D  --> p + 4He     3He(d,p)4He
    extrapolate: bool
        if True, T_ion outside of valid energy range (according to paper)
        will be extrapolated; if false those values will be set to NaN
    silent: bool
        if True, some (useful ?) output will be printed to console

    Returns
    -------
    float
        fusion reactivity in m^3/s
    '''

    func_name = 'get_fusion_reactivity_Hively'

    reaction_str = reaction_int2str( reaction, silent=True )

    # valid energy range according to paper
    energy_range = [1, 80]

    if not silent:
        print( '{0}:'.format(func_name) )
        print( '    fusion reactivity as obtained from the following paper:' )
        print( '    L.M. Hively, Nuclear Fusion, Vol. 17, No. 4 (1977)' )
        print( '    (more info in doc-string)' )
        print( '    reaction {0:d} ({1}), T_ion = {2} keV'.format(reaction, reaction_str, T_ion) )

    # set the coefficients of the fit equation
    if reaction_str == 'DT' or reaction_str == 'TD':             
        # Table I in Hively NF 1977 paper
        a1 = -21.377692
        a2 = -25.204054
        a3 = -7.1013427e-2
        a4 = 1.9375451e-4
        a5 = 4.9246592e-6
        a6 = -3.9836572e-8
        r  = .2935
    elif reaction_str == 'DD_a':           
        # Table III in Hively NF1977 paper
        a1 = -15.511891
        a2 = -35.318711
        a3 = -1.2904737e-2
        a4 = 2.6797766e-4
        a5 = -2.9198685e-6
        a6 = 1.2748415e-8
        r  = .3735
    elif reaction_str == 'DD_b':
        # Tabel IV in Hively NF1977 paper
        a1 = -15.993842         
        a2 = -35.017640
        a3 = -1.3689787e-2
        a4 = 2.7089621e-4
        a5 = -2.9441547e-6
        a6 = 1.2841202e-8
        r  = .3725
    elif reaction_str == '3HeD' or reaction_str == 'D3He':
        # Table II in Hively NF1977 paper
        a1 = -27.764468
        a2 = -31.023898
        a3 = 2.7809999e-2
        a4 = -5.5321633e-4
        a5 = 3.0293927e-6
        a6 = -2.5233325e-9
        r  = .3597

    # Eq. (5), referred to as S_5 in the paper
    sigma_v = 1e-6 * np.exp( a1/T_ion**r + a2 + a3*T_ion + a4*T_ion**2 + a5*T_ion**3 + a6*T_ion**4 )

    if not silent:
        print( '    ==> <sigma*v> = {0} m^3/s'.format(sigma_v) )

    # check if T_ion is within valid energy range
    if (np.amin(T_ion) < energy_range[0]) or (np.amax(T_ion) > energy_range[1]):
        print( '{0}:'.format(func_name) )
        print( '    WARNING: T_ion is outside of valid energy range' )
        if extrapolate:
            print( '{0}extrapolate was set to True (data will be extrapolated)'.format( 
                   ' '*13) )
        else:
            print( '{0}extrapolate was set to False (data will be set to NaN)'.format( 
                   ' '*13) )

            # if T_ion is array, set all values outside of range to NaN
            if np.size(T_ion) > 1:
                sigma_v[ T_ion < energy_range[0] ] = np.nan
                sigma_v[ T_ion > energy_range[1] ] = np.nan
            else:
                sigma_v = np.nan

    return sigma_v
#;}}}


def get_fusion_reactivity_Bosch( T_ion, reaction=1, extrapolate=True, silent=True ):
#;{{{

    '''
    Return the fusion reactivity for various fusion reactions. 

    Tabulated values from the following reference are used:
        H.-S. Bosch and G.M. Hale, Nuclear Fusion, Vol. 32, No. 4 (1992)
        https://doi.org/10.1088/0029-5515/32/4/I07

    Parameters
    ----------
    T_ion: float
        ion temperature in keV, valid range 0.2-100 keV
    reaction: int
        defines reaction considered, default value is 1
        possible values are:
        1: T + D    --> n + 4He     T(d,n)4He
        2: D + D    --> p + T       D(d,p)T
        3: D + D    --> n + 3He     D(d,n)3He
        4: 3He + D  --> p + 4He     3He(d,p)4He
    extrapolate: bool
        if True, T_ion outside of valid energy range (according to paper)
        will be extrapolated; if false those values will be set to NaN
    silent: bool
        if True, some (useful ?) output will be printed to console

    Returns
    -------
    float
        fusion reactivity in m^3/s
    '''

    func_name = 'get_fusion_reactivity_Bosch'

    reaction_str = reaction_int2str( reaction, silent=silent )

    # valid energy range according to paper
    energy_range = [.2, 100]

    if not silent:
        print( '{0}:'.format(func_name) )
        print( '    fusion reactivity as obtained from the following paper:' )
        print( '    H.-S. Bosch and G.M. Hale, Nuclear Fusion, Vol. 32, No. 4 (1992)')
        print( '    (more info in doc-string)' )
        print( '    reaction {0:d} ({1}), T_ion = {2} keV'.format(reaction, reaction_str, T_ion) )

    if reaction_str == 'DT' or reaction_str == 'TD':
        b_G     = 34.3827
        mr_c2   = 1124656.
        c1      = 1.17302e-9
        c2      = 1.51361e-2
        c3      = 7.51886e-2
        c4      = 4.60643e-3
        c5      = 1.35000e-2
        c6      = -1.06750e-4
        c7      = 1.36600e-5
        c7      = .0
    elif reaction_str == 'DD_a':
        b_G     = 31.3970
        mr_c2   = 937814.
        c1      = 5.65718e-12
        c2      = 3.41267e-3
        c3      = 1.99167e-3
        c4      = .0
        c5      = 1.05060e-5
        c6      = .0
        c7      = .0
    elif reaction_str == 'DD_b':
        b_G     = 31.3970
        mr_c2   = 937814.
        c1      = 5.43360e-12
        c2      = 5.85778e-3
        c3      = 7.68222e-3
        c4      = .0
        c5      = -2.96400e-6
        c6      = .0
        c7      = .0
    elif reaction_str == '3HeD' or reaction_str == 'D3He':
        b_G     = 68.7508
        mr_c2   = 1124572.
        c1      = 5.51036e-10
        c2      = 6.41918e-3
        c3      = -2.02896e-3
        c4      = -1.91080e-5
        c5      = 1.35776e-4
        c6      = .0
        c7      = .0

    # Eq. (13)
    theta = T_ion / ( 1. - T_ion*(c2+T_ion*(c4+T_ion*c6)) / (1.+T_ion*(c3+T_ion*(c5+T_ion*c7))) )

    # Eq. (14)
    chi   = ( b_G**2/(4.*theta) )**(1./3.)

    # reactivity as given in the paper in units of cm^3/s (with T_ion in keV)
    # Eq. (12)
    sigma_v = c1 * theta * np.sqrt( chi/(mr_c2*T_ion**3) ) * np.exp(-3.*chi)

    # scale to m^3/s
    sigma_v *= 1e-6

    if not silent:
        print( '    ==> <sigma*v> = {0} m^3/s'.format(sigma_v) )

    # check if T_ion is within valid energy range
    if (np.amin(T_ion) < energy_range[0]) or (np.amax(T_ion) > energy_range[1]):
        print( '{0}:'.format(func_name) )
        print( '    WARNING: T_ion is outside of valid energy range' )
        if extrapolate:
            print( '{0}extrapolate was set to True (data will be extrapolated)'.format( 
                   ' '*13) )
        else:
            print( '{0}extrapolate was set to False (data will be set to NaN)'.format( 
                   ' '*13) )

            # if T_ion is array, set all values outside of range to NaN
            if np.size(T_ion) > 1:
                sigma_v[ T_ion < energy_range[0] ] = np.nan
                sigma_v[ T_ion > energy_range[1] ] = np.nan
            else:
                sigma_v = np.nan

    return sigma_v
#;}}}


def log_interp1d(xVals, yVals, kind='linear'):
#;{{{
    """
    Wrapper around scipy.interpolate.interp1 for logarithmic datasets.

    Parameters
    ----------
    xVals: np.array
    yVals: np.array

    Returns
    -------
    log_interp: function 
    """

    logx = np.log10(xVals)
    logy = np.log10(yVals)

    # interpolate a 1D function to log10 of input-data
    f_interp_lin = interp.interp1d( logx, logy, kind=kind )

    # transform log back to linear scale
    f_interp_log = lambda x_new: np.power(10., f_interp_lin(np.log10(x_new)) )

    return f_interp_log
#;}}}


def get_fusion_reactivity_McNally( T_ion, reaction=1, extrapolate=True, silent=True ):
#;{{{
    '''
    Return the fusion reactivity for various fusion reactions.

    Fusion reactivity from tabular values (including interpolations to experimental values)
    Ref: J. Rand McNally, Fusion Reactivity Graphs and Tables for
         Charged Particle Reactions, ORNL/TM-6914, 1979
         https://doi.org/10.2172/5992170
    
    Parameters
    ----------
    T_ion: float
        ion temperature in keV, valid range 1-1000 keV
    reaction: int
        defines reaction considered, default value is 1
    extrapolate: bool
        if True, T_ion outside of valid energy range (according to paper)
        will be extrapolated; if false those values will be set to NaN
    silent: bool
        if True, some (useful ?) output will be printed to console

    Returns
    -------
    float
        fusion reactivity in m^3/s
    '''

    func_name = 'get_fusion_reactivity_McNally'

    reaction_str = reaction_int2str( reaction, silent=silent )

    if not silent:
        print( '{0}:'.format(func_name) )
        print( '    fusion reactivity as obtained from the following paper:' )
        print( '    J. Rand McNally, ORNL/TM-6914 (1979)')
        print( '    (more info in doc-string)' )
        print( '    reaction {0:d} ({1}), T_ion = {2} keV'.format(reaction, reaction_str, T_ion) )

    T_ion_tabulated = np.array( [  1, 2, 3, 4, 5, 6, 7, 8, 9
                                 ,10, 20, 30, 40, 50, 60, 70, 80, 90
                                 ,100, 200, 300, 400, 500, 600, 700, 800, 900
                                 ,1000
                                ] )

    energy_range = np.array( [ np.amin(T_ion_tabulated), np.amax(T_ion_tabulated) ] )

    if reaction_str == 'DD_b':       
        # Table I (page 9, top)
        sigma_v = np.array( [
              9.65e-29, 3.04e-27, 1.57e-26, 4.37e-26, 8.97e-26, 1.55e-25, 2.39e-25, 3.42e-25, 4.62e-25
            , 5.99e-25, 2.65e-24, 5.44e-24, 8.54e-24, 1.10e-23, 1.50e-23, 1.82e-23, 2.13e-23, 2.44e-23
            , 2.74e-23, 5.32e-23, 7.33e-23, 8.96e-23, 1.03e-22, 1.15e-22, 1.25e-22, 1.34e-22, 1.42e-22
            , 1.48e-22 
            ] )
    elif reaction_str == 'DD_a':     
        # Table I (page 9, bottom)
        sigma_v = np.array( [
              9.66e-29, 3.05e-27, 1.57e-26, 4.35e-26, 8.90e-26, 1.53e-25, 2.35e-25, 3.33e-25, 4.48e-25
            , 5.76e-25, 2.41e-24, 4.76e-24, 7.28e-24, 9.84e-24, 1.24e-23, 1.49e-23, 1.73e-23, 1.97e-23
            , 2.21e-23, 4.29e-23, 6.00e-23, 7.45e-23, 8.70e-23, 9.75e-23, 1.06e-22, 1.13e-22, 1.18e-22
            , 1.22e-22
            ] )
    elif reaction_str == 'DT' or reaction_str == 'TD':     
        # Table I (page 10, top)
        sigma_v = np.array( [
              6.27e-27, 2.83e-25, 1.81e-24, 5.86e-24, 1.35e-23, 2.53e-23, 4.14e-23, 6.17e-23, 8.57e-23
            , 1.13e-22, 4.31e-22, 6.65e-22, 7.93e-22, 8.54e-22, 8.76e-22, 8.76e-22, 8.64e-22, 8.46e-22
            , 8.24e-22, 6.16e-22, 4.90e-22, 4.13e-22, 3.63e-22, 3.28e-22, 3.02e-22, 2.83e-22, 2.68e-22
            , 2.55e-22
            ] )
    elif reaction_str == 'TT':     
        # Table I (page 10, bottom)
        sigma_v = np.array( [
              3.28e-29, 1.68e-27, 1.07e-26, 3.35e-26, 7.42e-26, 1.35e-25, 2.14e-25, 3.12e-25, 4.27e-25
            , 5.57e-25, 2.37e-24, 4.59e-24, 6.87e-24, 9.12e-24, 1.13e-23, 1.34e-23, 1.55e-23, 1.75e-23
            , 1.95e-23, 4.267e-23, 6.337e-23, 7.679e-23, 8.384e-23, 8.655e-23, 8.654e-23, 8.454e-23, 8.242e-23
            , 7.943e-23
            ] )
    elif reaction_str == 'T3He' or reaction_str == '3HeT':
        # Table I (page 11, top)
        sigma_v = np.array( [
              .0, .0, .0, .0, .0, .0, .0, .0, .0
            , 1.156e-26, 2.623e-25, 1.134e-24, 2.805e-24, 5.287e-24, 8.515e-24, 1.24e-23, 1.685e-23, 2.178e-23
            , 2.712e-23, 9.177e-23, 1.606e-22, 2.256e-22, 2.852e-22, 3.397e-22, 3.895e-22, 4.352e-22, 4.772e-22
            , 5.161e-22
            ] )
    elif reaction_str == 'D3He' or reaction_str == '3HeD':
        # Table I (page 11, bottom)
        sigma_v = np.array( [
              3.10e-32, 1.41e-29, 2.73e-28, 1.75e-27, 6.46e-27, 1.73e-26, 3.76e-26, 7.15e-26, 1.23e-25
            , 1.97e-25, 3.26e-24, 1.32e-23, 3.09e-23, 5.37e-23, 7.88e-23, 1.04e-22, 1.27e-22, 1.48e-22
            , 1.67e-22, 2.52e-22, 2.60e-22, 2.52e-22, 2.42e-22, 2.33e-22, 2.24e-22, 2.18e-22, 2.12e-22
            , 2.07e-22
            ] )
    elif reaction_str == '3He3He':
        sigma_v = np.array( [
              6.248073e-42, 3.553074e-37, 7.679157e-35, 1.434234e-33, 1.212912e-32, 6.463952e-32, 2.623214e-31
            , 8.195594e-31, 2.138978e-30, 4.860153e-30, 5.448750e-28, 4.924700e-27, 1.963613e-26, 5.186932e-26
            , 1.073749e-25, 1.900601e-25, 3.046649e-25, 4.511682e-25, 6.307594e-25, 4.073524e-24, 9.918522e-24
            , 1.774833e-23, 2.783765e-23, 4.056615e-23, 5.662810e-23, 7.607605e-23, 9.090703e-23, 1.251186e-22
            ] )
    elif reaction_str == 'p11B' or reaction_str == '11Bp':
        # Table I (page 16, top)
        sigma_v = np.array( [
              2.22553e-48, 1.23679e-37, 4.84893e-34, 3.48445e-32, 5.14265e-31, 3.43360e-30, 1.45314e-29
            , 4.64382e-29, 1.24057e-28, 2.93733e-28, 4.71713e-26, 3.96257e-25, 1.40318e-24, 3.50303e-24
            , 7.15115e-24, 1.26619e-23, 2.01130e-23, 2.93569e-23, 4.00945e-23, 1.62612e-22, 2.39482e-22
            , 2.79706e-22, 3.03366e-22, 3.19657e-22, 3.32697e-22, 3.44191e-22, 3.54833e-22, 3.64859e-22
            ] )

    # perform PCHIP 1D monotonic cubic interpolation
    #f_interp_sigmav = interp.PchipInterpolator( T_ion_tabulated, sigma_v )
    f_interp_sigmav = log_interp1d(T_ion_tabulated, sigma_v, kind='linear')
    sigma_v         = f_interp_sigmav( T_ion )

    # check if T_ion is within valid energy range
    if (np.amin(T_ion) < energy_range[0]) or (np.amax(T_ion) > energy_range[1]):
        print( '{0}:'.format(func_name) )
        print( '    WARNING: T_ion is outside of valid energy range' )
        if extrapolate:
            print( '{0}extrapolate was set to True (data will be extrapolated)'.format( 
                   ' '*13) )
        else:
            print( '{0}extrapolate was set to False (data will be set to NaN)'.format( 
                   ' '*13) )

            # if T_ion is array, set all values outside of range to NaN
            if np.size(T_ion) > 1:
                sigma_v[ T_ion < energy_range[0] ] = np.nan
                sigma_v[ T_ion > energy_range[1] ] = np.nan
            else:
                sigma_v = np.nan

    return sigma_v

    #return np.array( [ T_ion_tabulated, sigma_v ] )

#;}}}


def get_fusion_reactivity_Angulo( T_ion, reaction=9, reaction_str='', extrapolate=True, silent=True ):
#;{{{
    """
    Return the fusion reactivity for various fusion reactions.

    Fusion reactivity from analytical approximations, huge variety of reactions.
    Ref: C. Angulo et al., Nuclear Physics A 656 (1999) 3-183
         https://doi.org/10.1016/S0375-9474(99)00030-5
    
    Parameters
    ----------
    T_ion: float
        ion temperature in keV, valid range approximately .1-1000 keV
    reaction: int
        defines reaction considered, default value is 9
    extrapolate: bool
        if True, T_ion outside of valid energy range (according to paper)
        will be extrapolated; if false those values will be set to NaN
    silent: bool
        if True, some (useful ?) output will be printed to console

    Returns
    -------
    float
        fusion reactivity in m^3/s
    """

    func_name = 'get_fusion_reactivity_Angulo'

    if len(reaction_str) == 0:
        reaction_str = reaction_int2str( reaction, silent=silent )

    # temperatures in Angulo's formula are in 10^9 K
    eV2K    = consts.e/consts.Boltzmann     
    T9      = T_ion*1e3 * eV2K * 1e-9

    # valid temperature is given in units of T9 in paper, tranform to keV
    energy_range = np.array( [ .001, 10. ] )*1e9 / eV2K * 1e-3

    if not silent:
        print( '{0}:'.format(func_name) )
        print( '    fusion reactivity as obtained from the following paper:' )
        print( '    C. Angulo et al., Nuclear Physics A 656 (1979) 3-183')
        print( '    (more info in doc-string)' )
        print( '    reaction {0:d} ({1}), T_ion = {2} keV, T9 = {3} K'.format(reaction, reaction_str, T_ion, T9) )

    if reaction_str == 'pp':
        A = np.array( [ 4.08e-15, -3.381, 3.82, 1.51, 0.144, -1.14e-2 ] )
    else:
        print( '{0}: ERROR, no such reaction'.format( func_name ) )
        return np.nan

    sigma_v = 1./consts.Avogadro * A[0]*T9**(-2./3.) * np.exp(A[1]*T9**(-1./3.)) * (
                1. + A[2]*T9 + A[3]*T9**2 + A[4]*T9**3 + A[5]*T9**4 )

    # scale to m^3/s
    sigma_v *= 1e-6

    if not silent:
        print( '    ==> <sigma*v> = {0} m^3/s'.format(sigma_v) )

    # check if T_ion is within valid energy range
    if (np.amin(T_ion) < energy_range[0]) or (np.amax(T_ion) > energy_range[1]):
        print( '{0}:'.format(func_name) )
        print( '    WARNING: T_ion is outside of valid energy range' )
        if extrapolate:
            print( '{0}extrapolate was set to True (data will be extrapolated)'.format( 
                   ' '*13) )
        else:
            print( '{0}extrapolate was set to False (data will be set to NaN)'.format( 
                   ' '*13) )

            # if T_ion is array, set all values outside of range to NaN
            if np.size(T_ion) > 1:
                sigma_v[ T_ion < energy_range[0] ] = np.nan
                sigma_v[ T_ion > energy_range[1] ] = np.nan
            else:
                sigma_v = np.nan

    return sigma_v

#;}}}


def get_fusion_reactivity_Atzeni( T_ion, reaction=9, reaction_str='', extrapolate=True, silent=True ):
#;{{{
    """
    Return the fusion reactivity for various fusion reactions.

    Fusion reactivity from analytical approximations (very nice book). 
    Ref: S. Atzeni, J. Meyer-ter-Vehn, "Intertial Fusion: Beam Plasma Interactions,
         Hydrodynamics, Dense Plasma Physics" (Oxford Science Publications, 2004)
         ISBN-10: 0199568014
         ISBN-13: 978-0199568017
    
    Parameters
    ----------
    T_ion: float
        ion temperature in keV, valid range a bit unclear, for pp the 
        values for Angulo-paper are taken as it is references: ~ .1-1000 keV
    reaction: int
        defines reaction considered, default value is 9
    extrapolate: bool
        if True, T_ion outside of valid energy range (according to paper)
        will be extrapolated; if false those values will be set to NaN
    silent: bool
        if True, some (useful ?) output will be printed to console

    Returns
    -------
    float
        fusion reactivity in m^3/s
    """

    func_name = 'get_fusion_reactivity_Atzeni'

    energy_range = np.array( [.1, 1000] )

    if len(reaction_str) == 0:
        reaction_str = reaction_int2str( reaction, silent=silent )

    if not silent:
        print( '{0}:'.format(func_name) )
        print( '    fusion reactivity as obtained from the following paper:' )
        print( '    Atzeni & Meyer-ter-Vehn: Physics of Intertial Fusion (2004)')
        print( '    (more info in doc-string)' )
        print( '    reaction {0:d} ({1}), T_ion = {2} keV'.format(reaction, reaction_str, T_ion) )

    if reaction_str == 'pp':
        A = np.array( [ 1.56e-37, -14.94, 0.044, 2.03e-4, 5e-7 ] )
    else:
        print( '{0}: ERROR, no such reaction'.format( func_name ) )
        return np.nan

    sigma_v = A[0]*T_ion**(-2./3.) * np.exp(A[1]*T_ion**(-1./3.)) * (
                1. + A[2]*T_ion + A[3]*T_ion**2 + A[4]*T_ion**2 )

    # scale to m^3/s
    sigma_v *= 1e-6

    # check if T_ion is within valid energy range
    if (np.amin(T_ion) < energy_range[0]) or (np.amax(T_ion) > energy_range[1]):
        print( '{0}:'.format(func_name) )
        print( '    WARNING: T_ion is outside of valid energy range' )
        if extrapolate:
            print( '{0}extrapolate was set to True (data will be extrapolated)'.format( 
                   ' '*13) )
        else:
            print( '{0}extrapolate was set to False (data will be set to NaN)'.format( 
                   ' '*13) )

            # if T_ion is array, set all values outside of range to NaN
            if np.size(T_ion) > 1:
                sigma_v[ T_ion < energy_range[0] ] = np.nan
                sigma_v[ T_ion > energy_range[1] ] = np.nan
            else:
                sigma_v = np.nan

    return sigma_v

#;}}}


def main():
#;{{{

    print

    # possible values for the languare are 'en' and 'de'
    lang    = 'en'

    if lang == 'de':
        xlabel  = r'Ionentemperatur $T$ in keV'
        ylabel  = r'Fusionsreaktivität $\langle\sigma v\rangle$ in (m$^3$/s)'
    else:
        xlabel  = r'ion temperature $T$ in keV'
        ylabel  = r'fusion reactivity $\langle\sigma v\rangle$ in (m$^3$/s)'

    # set-up plot
    # (width, heigth) in inches
    fig1 = plt.figure( figsize=(8,6) )
    ax1  = fig1.add_subplot( 1,1,1 )
   
    # if empty, plot will be put out to X-window
    plot_fname = 'fusion_reactivity_{0}.png'.format(lang)
    #plot_fname  = ''

    # set fusion reaction, see function reaction_int2str
    reaction = 1

    # T_ion in keV
    T_ion = np.linspace(1,1000,1000)

    plot_Hively     = False
    plot_Bosch      = False
    plot_McNally    = True
    plot_Angulo     = True

    write_datasource2plot   = False
    write_plotcredit        = True

    if plot_Hively:
        txt_ref_str = 'Hively fit'
        ax1.plot( T_ion, get_fusion_reactivity_Hively(T_ion, reaction=1), 
                  label='T(d,n)4He', linewidth=2 )
        ax1.plot( T_ion, get_fusion_reactivity_Hively(T_ion, reaction=2), 
                  label='D(d,p)T', linewidth=2 )
        ax1.plot( T_ion, get_fusion_reactivity_Hively(T_ion, reaction=3),  
                  label='D(d,n)3He', linewidth=2 )
        ax1.plot( T_ion, get_fusion_reactivity_Hively(T_ion, reaction=4),  
                  label='3He(d,p)4He', linewidth=2 )

    if plot_Bosch:
        txt_ref_str = 'Bosch fit'
        ax1.plot( T_ion, get_fusion_reactivity_Bosch(T_ion, reaction=1),
                  label='T(d,n)4He', linewidth=3 )
        ax1.plot( T_ion, get_fusion_reactivity_Bosch(T_ion, reaction=2),
                  label='D(d,p)T', linewidth=3 )
        ax1.plot( T_ion, get_fusion_reactivity_Bosch(T_ion, reaction=3),
                  label='D(d,n)3He', linewidth=3 )
        ax1.plot( T_ion, get_fusion_reactivity_Bosch(T_ion, reaction=4),
                  label='3He(d,p)4He', linewidth=3 )

    if plot_McNally:
        lw_McNally  = 3
        txt_ref_str = 'McNally dataset'

        ax1.plot( T_ion, get_fusion_reactivity_McNally(T_ion, reaction=1), 
                  label='D+T', linewidth=lw_McNally )
        ax1.plot( T_ion, ( get_fusion_reactivity_McNally(T_ion, reaction=2)
                          +get_fusion_reactivity_McNally(T_ion, reaction=3)),
                  label='D+D', linewidth=lw_McNally )
        #ax1.plot( T_ion, get_fusion_reactivity_McNally(T_ion, reaction=5), 
        #          label='T+T', linewidth=lw_McNally )
        ax1.plot( T_ion, get_fusion_reactivity_McNally(T_ion, reaction=4), 
                  label=r'D+$^3$He', linewidth=lw_McNally )
        ax1.plot( T_ion, get_fusion_reactivity_McNally(T_ion, reaction=7), 
                  label=r'$^3$He+$^3$He', linewidth=lw_McNally )
        #ax1.plot( T_ion, get_fusion_reactivity_McNally(T_ion, reaction=8), 
        #          label=r'p+$^{11}$B', linewidth=lw_McNally )

    if plot_Angulo:
        lw_Angulo   = 3
        txt_ref_str = 'Angelo fit'

        ax1.plot( T_ion, get_fusion_reactivity_Angulo(T_ion, reaction=9),
                  label='p+p', linewidth=lw_Angulo )
        

    if write_datasource2plot:
        txt_ref_x0 = .75
        txt_ref_y0 = .84
        fig1.text( txt_ref_x0, txt_ref_y0, txt_ref_str, fontsize=10 )

    ax1.set_xlim( np.amin(T_ion), np.amax(T_ion) )
    if plot_Angulo:
        ax1.set_ylim( [1e-50, 1e-20] )
    else:
        ax1.set_ylim( [1e-26, 1e-20] )
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax1.set_xlabel( xlabel )
    ax1.set_ylabel( ylabel )

    ax1.grid()

    legend = ax1.legend( loc='best', fontsize=15 )
    legend.get_frame().set_alpha(0.5)

    # per default, write credits into plot, to ensure that people know they can use the plot
    # every plot appearing somewhere on the internet should contain information on how to 
    # use it, otherwise it is useless in terms of re-usability
    # you probably want to remove it when you make your own plot
    # (attribution would still be gratefully acknowledged :)
    # also note that the license refers only to that specific plot
    # the license for the code is mentioned above and in the LICENSE file
    if write_plotcredit:
        credit_str = u'{0}, CC BY-SA 4.0'.format( __author__ )
        fig1.text( .7, .885, credit_str, fontsize=7 )

    make_plot( plot_fname )

#;}}}


def test():
#;{{{

    print( 'testing' )

    reaction_int = 1
    reaction_str = reaction_int2str( reaction_int, silent=False )
    print( 'reaction = {0}, {1}'.format( reaction_int, reaction_str ) )

    T_ion   = 1.

    sigma_v_Hively = get_fusion_reactivity_Hively( T_ion, reaction=reaction_int, silent=False, extrapolate=False )
    sigma_v_Bosch  = get_fusion_reactivity_Bosch( T_ion, reaction=reaction_int, silent=False, extrapolate=True )
    sigma_v_McNally= get_fusion_reactivity_McNally( T_ion, reaction=1, silent=False )

    print( 'T_ion = {0} keV, sigma_v_Hively = {1:e}, sigma_v_Bosch = {2:e}, sigma_v_McNally = {3:e}'.format(
            T_ion, sigma_v_Hively, sigma_v_Bosch, sigma_v_McNally) )

    sigma_v_Angulo  = get_fusion_reactivity_Angulo( T_ion, reaction=9, reaction_str='', extrapolate=True, silent=False )
    sigma_v_Atzeni  = get_fusion_reactivity_Atzeni( T_ion, reaction=9, reaction_str='', extrapolate=True, silent=False )
    print( 'T_ion = {0} keV, sigma_v_Angulo = {1:e}, sigma_v_Atzeni = {2:e}'.format(
            T_ion, sigma_v_Angulo, sigma_v_Atzeni) )
#;}}}


if __name__ == '__main__':
#    test()
    main()

