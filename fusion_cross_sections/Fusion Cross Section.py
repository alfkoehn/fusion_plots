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


#%%
"""
Plot the total cross section in m^2 for various species vs incident energy in keV
"""
E = np.logspace(0, 3, 501)
barns_to_SI = 1e-24 * 1e-4 # in m^2

sigma_DD = barns_to_SI*(cross_section_NRL(E, 'DD_a') + cross_section_NRL(E, 'DD_b'))
sigma_DT = barns_to_SI*cross_section_NRL(E, 'DT')
sigma_DHe3 = barns_to_SI*cross_section_NRL(E, 'DHe3')

fig, ax = plt.subplots()
ax.loglog(E, sigma_DD, E, sigma_DT, E, sigma_DHe3, lw=3)

ax.set_ylim([1e-32, 2e-27])
ax.set_xlim(1,1e3)

# Format tick values as numbers
maj_formatter = matplotlib.ticker.ScalarFormatter()
maj_formatter.set_scientific(False)
ax.get_xaxis().set_major_formatter(maj_formatter) 
ax.set_xticks([1, 10, 100, 1000])

# Top x-axis for energies in millions of K
def keV_to_K(keV):
    return keV/k_B*e/1e3
def K_to_keV(K):
    return K*k_B/e*1e3

ax2 = ax.secondary_xaxis('top', functions=(keV_to_K, K_to_keV))
ax2.set_xticks([20, 100, 1000, 5000])
# force matplotlib LogFormatterSciNotation to not use scientific notation until 10^5
plt.rcParams['axes.formatter.min_exponent'] = 5

[a.tick_params(labelsize=14) for a in (ax, ax2)]

ax.grid(True, which='minor')
ax.set_xlabel('Deuteron Energy [keV]', fontsize=16)
ax.set_ylabel('Cross section $\sigma$ [$m^2$]', fontsize=16)
ax.legend(('D-D', 'D-T', 'D-He$^3$'), loc='best', fontsize=18)
ax2.set_xlabel('$T$ [million K]', fontsize=16)

fig.tight_layout()
# fig.text( .71, .98, credit_str, fontsize=7)
fig.savefig('cross_sections_vs_temperature.png', dpi=300)

#%%
"""
Plot the Reaction rates (reactivity) <sigma v> in cm^3/s as a function of E,
the energy in keV of the incident particle [the first ion of the reaction label]
Data from NRL formulary.
See also
https://en.wikipedia.org/wiki/Nuclear_fusion#Maxwell_averaged_nuclear_cross_sections
"""
E, DT, DD, DHe3,TT,THe3 = np.loadtxt('reaction_rates_vs_energy_incident_particle.csv', 
                   skiprows=1, unpack=True)

fig, ax = plt.subplots()
ax.loglog(E, DD, 
          E, DT, 
          E, DHe3, lw=3)
ax.grid(True)
ax.set_xlabel('Temperature [keV]', fontsize=16)
ax.set_ylabel(r'Reaction Rates $\left<\sigma v\right>$ [$cm^3/s$]', fontsize=16)
ax.tick_params(labelsize=14)
ax.set_ylim([1e-26, 2e-15])
ax.legend(('D-D', 'D-T', 'D-He$^3$'), loc='best', fontsize=18)
fig.tight_layout()
# fig.text( .71, .98, credit_str, fontsize=7)
fig.savefig('reaction_rates_section_vs_E.png', dpi=300)

