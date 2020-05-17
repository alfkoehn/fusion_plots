# -*- coding: utf-8 -*-
""" 
Plot the reaction rates and the cross section as a function of energy

Data taken from NRL Formulary:
https://www.nrl.navy.mil/ppd/content/nrl-plasma-formulary
See also The Magnetic Fusion Energy Formulary:
https://github.com/MFEFormulary/MFEFormulary
"""
__author__      = 'Julien Hillairet'
__email__       = 'julien.hillairet@cea.fr'
__copyright__   = 'CEA/IRFM'
__license__     = 'MIT'

# write credits into plot, to ensure that people know they can use the plot
# (somebody once told me, every plot appearing somewhere in the internet
#  should contain information on how to use it, otherwise it is useless)
# you probably want to remove it when you make you own plot
# attribution would still be gratefully acknowledged :)
# also note that the licence refers only to that specific plot
# the licence for the code is mentioned above and in the LICENCE file
credit_str = f'{__author__}, CC BY-SA 4.0'

# import standard modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
from scipy.constants import e, Boltzmann as k_B

#%%
def cross_section_NRL(E, reaction='DT'):
    """
    The total cross section in barns (1 barns=1e-24 cm^2) as a function of E, 
    the energy in keV of the incident particle.
    
    Formula from NRL Formulary
    https://www.nrl.navy.mil/ppd/content/nrl-plasma-formulary
    See also
    https://en.wikipedia.org/wiki/Nuclear_fusion#Formulas_of_fusion_cross_sections

    
    Parameters
    ----------
    E : array
        energy (in keV) of the incident particle towards a target ion at rest 
    reaction : str
        Reaction: 'DT', 'DD_a', 'DD_b', 'DHe3', 'TT', 'THe3'. Default is 'DT'.

    Returns
    -------
    sigma_v: array
        Total cross section in barns    

    """
    cross_section_coefficients = {
        'DD_a': [46.097, 372, 4.36e-4, 1.220, 0],
        'DD_b': [47.88, 482, 3.08e-4, 1.177, 0], 
        'DT':   [45.95, 50200, 1.368e-2, 1.076, 409],
        'DHe3': [89.27, 25900, 3.98e-3, 1.297, 647],
        'TT':   [38.39, 448, 1.02e-3, 2.09, 0],
        'THe3': [123.1, 11250, 0, 0, 0],
        }
    A = cross_section_coefficients[reaction]
    sigma_T = (A[4]+((A[3]-A[2]*E)**2+1)**(-1) * A[1])/(E*(np.exp(A[0]/np.sqrt(E))-1))
    return(sigma_T)


def cross_section_Miley( T_ion, reaction='DT') :
#;{{{
    """
    """

    sigma_T = (T*(np.exp(A[0]/np.sqrt(T)) - 1.))**(-1) * ((A[1]/(1.+(A[2]*T-A[3])**2)) + A[4] )
#;}}}


def cross_section_Bosch( T_ion, reaction='DT' ):
#;{{{
    """
    Bosch Nuclear Fusion (1992)

    Energy refers to the energy available in the center-of-mass
    frame (CM). For particle A with mass m_A striking a stationary
    particle B, following holds: E_A = E*(m_A+m_B)/m_B
    """

    # T_ion: keV
    # cross-section in mb (millibarn, corresponding to 1e-31 m^2)
    # reactions as a function of the energy in the centre-of-mass (CM) frame

    T = T_ion

    if reaction == 'DT' or reaction == 'TD':
        # valid energy range in keV
        energy_range = np.array( [.5, 550] )
        B_G = 34.3827   # sqrt(keV)
        A   = [ 6.927e4, 7.454e8, 2.050e6, 5.2002e4, .0 ]
        B   = [ 6.38e1, -9.95e-1, 6.981e-5, 1.728e-4 ]
    elif reaction == 'He3D' or reaction == 'DHe3' :
        energy_range = np.array( [.3, 900] )
        B_G = 68.7508   # sqrt(keV)
        A   = [ 5.7501e6, 2.5226e3, 4.5566e1, .0, .0 ]
        B   = [ -3.1995e-3, -8.5530e-6, 5.9014e-8, .0 ] 
    elif reaction == 'DD':
        sigma_T_a = cross_section_Bosch( T_ion, reaction='DT_a' )
        sigma_T_b = cross_section_Bosch( T_ion, reaction='DT_b' )
        return sigma_T_a + sigma_T_b
    elif reaction == 'DD_a':
        energy_range = np.array( [ .5, 5000] )
        B_G = 31.3970   # sqrt(keV)
        A   = [ 5.5576e4, 2.1054e2, -3.2638e-2, 1.4987e-6, 1.8181e-10 ]
        B   = [ .0, .0, .0, .0 ]
    elif reaction == 'DD_b':
        energy_range = np.array( [ .5, 4900] )
        B_G = 31.3970   # sqrt(keV)
        A   = [ 5.3701e4, 3.3027e2, -1.2706e-1, 2.9327e-5, -2.5151e-9 ]
        B   = [ .0, .0, .0, .0 ]

    # Padé polynomial
    S_func = lambda T: ( 
            A[0] + T*(A[1] + T*(A[2] + T*(A[3] + T*A[4]))) / 
            (1.  + T*(B[0] + T*(B[1] + T*(B[2] + T*B[3]))))
            )

    # cross section in mbarn
    sigma_T = S_func(T) / (T*np.exp(B_G/np.sqrt(T)))

    # set values outside of valid energy range to NaN
    sigma_T[ T < energy_range[0] ] = np.nan
    sigma_T[ T > energy_range[1] ] = np.nan

    # return cross-section in barn
    return sigma_T*1e-3

#;}}}

#%%
# keV <-> K conversions for upper x-axis
def keV_to_K(keV):
    return keV/k_B*e/1e3
def K_to_keV(K):
    return K*k_B/e*1e3


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
        #fig.tight_layout()
        #fig.savefig( fname_plot, dpi=600 )
        print( 'written plot into file {0}'.format(fname_plot) )
    else:
        plt.show()

#;}}}


def main():
#;{{{
    """
    Plot the total cross section in m^2 for various species vs incident energy in keV
    """

    dataset = 'Bosch'

    # ion temperature in units of keV 
    # (actually energy, but in plasma physics we just call it temperature :-)
    T_ion = np.logspace(0, 3, 501)

    # 1 barn = 1e-24 cm^2 = 1e-28 m^2
    barns_to_SI = 1e-28

    # get cross section
    if dataset == 'NRL':
        sigma_DD    = barns_to_SI*(cross_section_NRL(T_ion, 'DD_a') + cross_section_NRL(T_ion, 'DD_b'))
        sigma_DT    = barns_to_SI*cross_section_NRL(T_ion, 'DT')
        sigma_DHe3  = barns_to_SI*cross_section_NRL(T_ion, 'DHe3')
    elif dataset == 'Bosch':
        sigma_DD    = barns_to_SI*(cross_section_Bosch(T_ion, 'DD_a') + cross_section_NRL(T_ion, 'DD_b'))
        sigma_DT    = barns_to_SI*cross_section_Bosch(T_ion, 'DT')
        sigma_DHe3  = barns_to_SI*cross_section_Bosch(T_ion, 'DHe3')

    fig, ax1 = plt.subplots()
    ax1.loglog(T_ion, sigma_DD, T_ion, sigma_DT, T_ion, sigma_DHe3, lw=3)

    ax1.set_ylim([1e-32, 2e-27])
    ax1.set_xlim(1,1e3)

    # Create an upper x-axis with temperature in millions K
    # and format the logscale ticks in natural numbers instead of scientific notation
    # Format tick values as numbers
    maj_formatter = matplotlib.ticker.ScalarFormatter()
    maj_formatter.set_scientific(False)
    ax1.get_xaxis().set_major_formatter(maj_formatter) 
    ax1.set_xticks([1, 10, 100, 1000])
    ax2 = ax1.secondary_xaxis('top', functions=(keV_to_K, K_to_keV))
    ax2.set_xticks([20, 100, 1000, 5000])
    # force matplotlib LogFormatterSciNotation to not use scientific notation until 10^5
    plt.rcParams['axes.formatter.min_exponent'] = 5
    [a.tick_params(labelsize=14) for a in (ax1, ax2)]

    ax1.grid(True, which='both')
    ax1.set_xlabel('Deuteron Energy [keV]', fontsize=16)
    ax1.set_ylabel('Cross section $\sigma$ [$\mathrm{m}^2$]', fontsize=16)
    ax1.legend(('D-D', 'D-T', 'D-He$^3$'), loc='best', fontsize=18)
    ax2.set_xlabel('$T$ [million K]', fontsize=16)

    # set ticks to point inwards (looks nicer, in my opinion)
    ax1.tick_params(axis='both', which='both', direction='in', top=False, right=True)
    ax2.tick_params(axis='both', which='both', direction='in', top=True,  right=False)

    #fig.tight_layout()
    # fig.text( .71, .98, credit_str, fontsize=7)
    #fig.savefig('cross_sections_vs_temperature__{0}.png'.format(dataset), dpi=300)
    make_plot( fname_plot='cross_sections_vs_temperature__2.png' )
#;}}}


def main_test():

    T_ion = np.array( [3, 4, 5, 6, 7, 8, 9, 10, 100, 400 ] )

    for ii in range(len(T_ion)):
        print( ' T_ion = {0:3.0f} keV  ==>  sigma = {1:9.3e}'.format(
               T_ion[ii], cross_section_Bosch( T_ion[ii], reaction='DT' ) ) )
#               T_ion[ii], cross_section_Bosch( T_ion[ii], reaction='He3D' ) ) )
#               T_ion[ii], cross_section_Bosch( T_ion[ii], reaction='DD_a' ) ) )
#               T_ion[ii], cross_section_Bosch( T_ion[ii], reaction='DD_b' ) ) )

if __name__ == '__main__':
    main()
    #main_test()

