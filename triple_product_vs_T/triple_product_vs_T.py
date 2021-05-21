# coding=utf-8

"""
Plot the triple product as a function of T_ion for various experiments.

Simply run this script to produce a png plot:
    $ python triple_product_vs_T.py
"""

__author__      = 'Alf Köhn-Seemann'
__email__       = 'koehn@igvp.uni-stuttgart.de'
__copyright__   = 'University of Stuttgart'
__license__     = 'MIT'

# import standard modules
import numpy as np
import matplotlib.pyplot as plt

# change some default properties of matplotlib
plt.rcParams.update({'font.size':12})
# force ticks to point inwards
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top']       = True
plt.rcParams['ytick.right']     = True


def get_DT_fusion_reactivity( T_ion, silent=True ):
#;{{{
    '''
    Return the fusion reactivity for the D-T fusion reaction. 

    Tabulated values from the following reference are used:
        H.-S. Bosch and G.M. Hale, Nuclear Fusion, Vol. 32, No. 4 (1992)
        https://doi.org/10.1088/0029-5515/32/4/I07
    The values were extracted manually using the very nice tool WebPlotDigitizer.

    Parameters
    ----------
    T_ion: float
        ion temperature in keV, valid range 0.2-100 keV
    silent: bool
        if True, some (useful ?) output will be printed to console

    Returns
    -------
    float
        fusion reactivity in m^3/s
    '''

    if not silent:
        print( 'get_DT_fusion_reactivity' )
        print( '    check the doc-string for the corresponding reference' )

    # the paper by Bosch contains data for more fusion reactions
    # here, we focus on DT-fusion only
    # reaction = 1 refers to T(d,n)4He
    reaction = 1

    if reaction == 1:
        b_G     = 34.3827
        mr_c2   = 1124656.
        c1      = 1.17302e-9
        c2      = 1.51361e-2
        c3      = 7.51886e-2
        c4      = 4.60643e-3
        c5      = 1.35000e-2
        c6      = -1.06750e-4
        c7      = 1.36600e-5

    theta = T_ion / ( 1. - T_ion*(c2+T_ion*(c4+T_ion*c6)) / (1.+T_ion*(c3+T_ion*(c5+T_ion*c7))) )
    chi   = ( b_G**2/(4.*theta) )**(1./3.)

    # reactivity as given in the paper in units of cm^3/s (with T_ion in keV)
    sigma_v = c1 * theta * np.sqrt( chi/(mr_c2*T_ion**3) ) * np.exp(-3.*chi)
    # scale to m^3/s
    sigma_v *= 1e-6

    if not silent:
        print( '    T_i = {0:5.1f} keV => (D+T) reactivity = {1:8.3e} m^3/s',format(T_ion, sigma_v) )

    return sigma_v
#;}}}


def get_experimental_dataset( dataset='Bosch', silent=True ):
#;{{{
    '''
    Return experimental values of triple-product achieved in different devices. 

    Different datasets can be returned, based on the following references:
        'Bosch': Plasma Physics: Confinement, Transport and Collective Effects (2005), 
                 Editor Dinklage, Chapter "Nuclear Fusion" H.-S. Bosch
                 https://doi.org/10.1007/11360360_17
                 Data manually extracted from the plot using the very
                 useful tool WebPlotDigitizer.
        'my_dset': separate datapoints extracted from various references,
                   more details are included in my (hand-written) notebook.
                   Following list shows the refs used for each experiment:
                   START: Gates et al., PoP Vol. 5, No. 5, p 1775 (1998)
                          https://doi.org/10.1063/1.872819
                   Globus-M: Gusev et al., NF Vol. 55, p 104016 (2015)
                             https://doi.org/10.1088/0029-5515/55/10/104016
                   NSTX: Kaye et al., NF Vol. 55, p 104002 (2015)
                         https://doi.org/10.1088/0029-5515/55/10/104002
                   MAST: Chapman et al., NF Vol. 55, p 104008 (2015)
                         https://doi.org/10.1088/0029-5515/55/10/104008
                   W7-X (lim): Hirsch et al., NF Vol. 52, p 086010 (2017)
                               https://doi.org/10.1088/1741-4326/aa7372
                   W7-X (div): Pedersen et al., PPCF Vol. 61, p 014035 (2019)
                               https://doi.org/10.1088/1361-6587/aaec25
                   ITER: ########## REF STILL MISSING #####################
                   TCV: Karpuskov et al., FED Vol. 123, p 468 (2017)
                        https://doi.org/10.1016/j.fusengdes.2017.02.076
        'EUROfusion': data manually extracted from plot used by EUROfusion
                      using the very useful tool WebPlotDigitizer.

    Parameters
    ----------
    dataset: str
        possible values are 'Bosch', 'my_dset', 'EUROfusion'
         
    silent: bool
        if True, some (useful ?) output will be printed to console

    Returns
    -------
    tuple (of numpy-arrays)
        T_vals: ion temperature in keV, 
        nTtau_vals: triple product in m^-3 keV s, 
        device_types: type of device (1:stellarator, 2:tokamak, 3:spherical tokamak), 
        names: name of device
    '''

    if not silent:
        print( 'get_experimental_dataset' )
        print( '    dataset chosen: {0}'.format(dataset) )
        print( '    check the doc-string for the corresponding references' )

    if dataset == 'Bosch':
        names = np.array( [ 'T3', 'T3', 'Pulsator', 'T10', 
                            'ASDEX', 'ASDEX', 'Tore Supra', 'ASDEX Upgrade', 
                            'JT60', 'Alcator C-mod', 'TFTR', 'Alcator', 
                            'JET', 'DIII-D', 'ASDEX Upgrade', 'DIII-D', 
                            'TFTR', 'TFTR (DT)', 'JT60U', 'JT60U', 
                            'JET (DT)', 'JET', 'JET (DT)', '?', 
                            'JT60U', 'LHD', 'LHD', 'W7-AS', 
                            'W7-A'
                          ] )
        device_types = np.array( [2, 2, 2, 2, 
                                  2, 2, 2, 2, 
                                  2, 2, 2, 2, 
                                  2, 2, 2, 2, 
                                  2, 2, 2, 2, 
                                  2, 2, 2, 0, 
                                  2, 1, 1, 1, 
                                  1
                                 ] )
        T_vals = np.array( [ 0.16974119, 0.30918382, 0.22998638, 0.73804885,
                             4.96777849, 2.2553371 , 2.57170267, 3.35038255,
                             6.73158083, 1.50030032, 1.49146486, 3.07362999,
                             4.42561292, 6.33446757, 11.0294331, 12.6267789,
                             34.1913846, 40.4518789, 16.190024,  15.939778,
                             20.5655601, 28.3716432, 33.7292375, 40.408576,
                             44.9754023, 1.50931078, 9.16780799, 0.83739708,
                             0.61829778
                           ] )
        nTtau_vals = np.array( [1.75839569e+16, 2.10183008e+17, 8.34943817e+17, 3.05070787e+18,
                                5.39415148e+18, 1.60864212e+19, 2.79836303e+19, 4.47870397e+19,
                                3.08640860e+19, 7.40596733e+19, 1.80870567e+20, 8.34919183e+19,
                                1.62570081e+20, 1.69728539e+20, 8.37008333e+19, 2.06716580e+20, 
                                2.71606706e+20, 5.15793728e+20, 6.63979351e+20, 8.80494689e+20,
                                7.44951684e+20, 8.82853215e+20, 1.16051029e+21, 1.24463067e+21,
                                1.54281713e+21, 1.55385139e+19, 2.11990093e+19, 1.66398614e+18,
                                5.83188960e+17
                               ] )
    elif dataset == 'my_dset':
        names = np.array( ['START', 'Globus-M', 'NSTX', 'MAST', 
                           'W7-X (lim)', 'W7-X (div)', 'ITER', 'TCV'
                          ] )
        device_types = np.array( [ 3, 3, 3, 3, 
                                   1, 1, 2, 2
                                 ] )
        T_vals = np.array( [ 0.25, 0.5, 1., 0.9, 
                             1.,   3.5, 20., 3.7
                           ] )
        nTtau_vals = np.array( [ 2.50000000e+16, 4.50000000e+16, 6.00000000e+17, 8.10000000e+17,
                                 2.00000000e+18, 6.60000000e+19, 3.00000000e+21, 3e18
                               ] )
    elif dataset == 'EUROfusion':
        names = np.array( [ 'T3', 'TFR', 'T10', 'PLT', 
                            'ASDEX', 'ALC-A', 'TFR', 'ASDEX',
                            'TEXTOR', 'PLT', 'Tore Supra', 'FT', 
                            'ALC-C', 'TFTR', 'JT-60', 'DIII-D', 
                            'JET', 'ASDEX-U', 'DIII-D', 'TFTR', 
                            'JT-60U', 'JET', 'DIII-D', 'TFTR', 
                            'JET', 'JET', 'JET', 'JET', 
                            'TFTR', 'JT-60U', 'ITER'
                          ] )
        device_types = np.array( [ 2, 2, 2, 2, 
                                   2, 2, 2, 2, 
                                   2, 2, 2, 2, 
                                   2, 2, 2, 2, 
                                   2, 2, 2, 2, 
                                   2, 2, 2, 2, 
                                   2, 2, 2, 2, 
                                   2, 2, 2
                                 ] )
        T_vals = np.array( [ 0.29457498, 0.87620758, 0.77754643, 0.95833535,
                             0.67219729, 0.72430616, 1.6902491 , 2.3214469,
                             4.68336189, 6.95684699, 2.76667723, 0.99478727,
                             1.29670285, 1.52249236, 2.72567086, 5.04641705,
                             5.1414959 , 9.91820781, 15.6396948, 19.3482071,
                             9.88125103, 17.956236 , 21.1617331, 29.611856,
                             11.3872861, 18.294547 , 31.084376 , 35.0286038,
                             47.3964468, 47.7516439, 19.2761126
                           ] )
        nTtau_vals = np.array( [1.90726752e+17, 1.05386684e+18, 2.35462683e+18, 3.91227391e+18,
                                7.38011641e+18, 1.33450920e+19, 1.41714706e+18, 1.48969763e+19,
                                1.36882186e+19, 3.44589835e+18, 3.62231044e+19, 3.97568513e+19,
                                9.18858856e+19, 1.76292161e+20, 1.00849818e+20, 1.56596240e+20,
                                2.08802354e+20, 1.12577519e+20, 7.37387379e+19, 1.60622614e+20,
                                3.41108325e+20, 3.41108325e+20, 4.62590403e+20, 3.49878839e+20,
                                8.43587429e+20, 9.33751973e+20, 7.36763646e+20, 9.74108044e+20,
                                9.02674252e+20, 1.57793438e+21, 3.86947102e+21
                               ] )

    return T_vals, nTtau_vals, device_types, names

#;}}}


def make_plot( fname_plot='' ):
#;{{{
    '''
    Output a plot, either to X-window (default) or into file. 

    Parameters
    ----------
    fname_plot: str
        possible values are 'Bosch', 'my_dset', 'EUROfusion'

    Returns
    -------
    '''


    if len(fname_plot) > 0:
        plt.savefig( fname_plot, dpi=600, bbox_inches='tight' )
        print( 'written plot into file {0}'.format(fname_plot) )
    else:
        plt.show()

#;}}}


def main():
#;{{{

    print( '\nTo make a plot of n*T*tau_E as a function of T, data from various' )
    print( 'references is used. Check the doc-string of each function to get a' )
    print( 'list of the corresponding references.')
    print( 'Let me know if you have questions, requests or found some bugs.')
    print( '                                    -- Alf Köhn-Seemann, April 2020\n' )

    # set fname_plot to an empty string to plot into X-window
    fname_plot = 'triple_product_vs_T.png'

    # load Bosch dataset
    T_vals__B, nTtau_vals__B, device_types__B, names__B = get_experimental_dataset( dataset='Bosch' )
    # load my own dataset
    T_vals__my, nTtau_vals__my, device_types__my, names__my = get_experimental_dataset( dataset='my_dset' )
    # combine datasets
    T_vals          = np.append( T_vals__B,         T_vals__my )
    nTtau_vals      = np.append( nTtau_vals__B,     nTtau_vals__my )
    device_types    = np.append( device_types__B,   device_types__my )
    names           = np.append( names__B,          names__my )

    # constants required for small lambda functions below
    c_br    = 1.04e-19          # m^3 eV^1/2 s^-1
    E_alpha = 3.52e6            # MeV

    # ignition: alpha-heating used to sustain fusion reaction with no external heating
    F_ignition  = lambda T_ion, Z_eff: 12.*T_ion**2 / ( get_DT_fusion_reactivity(T_ion*1e-3)*E_alpha
                                                       - 4.*c_br*Z_eff*np.sqrt(T_ion) )

    # bremsstrahlung limit
    # requirement: bremsstrahlung losses < total energy loss rate per unit volume
    #   => P_rad <= W/tau_E
    # re-arranging yields n*T*tau_E <= 3*T^1.5/c_br
    F_bremslimit = lambda T_ion, Z_eff: 3.*T_ion**(1.5) / (Z_eff*c_br)

    # set-up plot
    fig1 = plt.figure( figsize=(8,6) )      # (width, height)
    ax1  = fig1.add_subplot( 1,1,1 )

    # filled area indicating bremsstrahlung limit
    T_full = np.logspace( np.log10(.1e3), np.log10(100e3), 100 )
    ax1.fill_between( x=T_full*1e-3, y1=F_bremslimit(T_full,1)*1e-3, y2=1e22, 
                      color='grey',
                    )
    ax1.annotate( 'bremsstrahlung limit', xy=( 0.15, 1.1e21), color='.85', rotation=29.5 )

    # filled area indicating ignition
    T_ignition  = np.logspace( np.log10(1e3), np.log10(100e3), 100 )
    F_ign_Z1    = F_ignition( T_ignition, 1 )
    ax1.fill_between( T_ignition[ F_ign_Z1>0 ]*1e-3, 
                      F_ign_Z1[ F_ign_Z1>0 ]*1e-3 , 
                      1e22 )
    ax1.annotate( 'ignition (DT)', xy=( 8.1, 5.3e21), color='.85' )

    # plot experimental values
    # stellarators
    ax1.plot( T_vals[ device_types == 1 ], nTtau_vals[ device_types == 1 ],
              marker='o', color='red', linestyle='None',
              label='Stellarator'
            )
    # tokamaks
    ax1.plot( T_vals[ device_types == 2 ], nTtau_vals[ device_types == 2 ],
              marker='s', color='blue', linestyle='None',
              label='Tokamak'
            )
    # spherical tokamaks
    ax1.plot( T_vals[ device_types == 3 ], nTtau_vals[ device_types == 3 ],
              marker='D', color='green', linestyle='None',
              label='Spherical tokamak'
            )

    # write annotations to experimental values
    # note: using some external graphics program (like inkscape, gimp)
    #       would certainly be faster than the following way...
    # list with cases which cannot be annotated automatically
    # because labels would otherwise overlap with something
    extra_labels = [ 'TFTR (DT)', 'JT60U', 'Alcator C-mod', 'Alcator', 'JET (DT)', 
                     '?', 'Tore Supra', 'LHD', 'JET', 'ASDEX Upgrade', 'DIII-D',
                     'ASDEX', 'W7-A', 'W7-AS' ]
    if len(fname_plot) > 0:
        dpi_scale   = 6
    else:
        dpi_scale   = 1
    for ii in range( len(names) ):
        # automatic annotation
        if names[ii] not in extra_labels:
            plt.annotate( names[ii], 
                          xy=( T_vals[ii], nTtau_vals[ii]), 
                          xytext=(8*dpi_scale,-5*dpi_scale), 
                          textcoords='offset pixels',
                        )
        # handles the above defined extra cases
        elif (names[ii] == 'Alcator C-mod') or (names[ii] == 'Tore Supra'):
            plt.annotate( names[ii], 
                          xy=( T_vals[ii], nTtau_vals[ii]), 
                          xytext=(-9.8*dpi_scale*len(names[ii]),-6*dpi_scale), 
                          textcoords='offset pixels',
                        )
        elif names[ii] == 'W7-A':
            plt.annotate( names[ii], 
                          xy=( T_vals[ii], nTtau_vals[ii]), 
                          xytext=(-8*dpi_scale,-20*dpi_scale), 
                          textcoords='offset pixels',
                        )
        elif names[ii] == 'W7-AS':
            plt.annotate( names[ii], 
                          xy=( T_vals[ii], nTtau_vals[ii]), 
                          xytext=(-12*dpi_scale*len(names[ii]),-5*dpi_scale), 
                          textcoords='offset pixels',
                        )
        elif names[ii] == 'Alcator':
            plt.annotate( names[ii], 
                          xy=( T_vals[ii], nTtau_vals[ii]), 
                          xytext=(4.*dpi_scale,+2*dpi_scale), 
                          textcoords='offset pixels',
                        )
        elif names[ii] == 'TFTR (DT)':
            plt.annotate( names[ii], 
                          xy=( T_vals[ii], nTtau_vals[ii]), 
                          xytext=(-6*dpi_scale,5*dpi_scale), 
                          textcoords='offset pixels',
                        )
        elif names[ii] == 'LHD':
            if T_vals[ii] < 2:
                plt.annotate( names[ii], 
                              xy=( T_vals[ii], nTtau_vals[ii]), 
                              xytext=(-14*dpi_scale*len(names[ii]),-6*dpi_scale), 
                              textcoords='offset pixels',
                            )
            else:
                plt.annotate( names[ii], 
                              xy=( T_vals[ii], nTtau_vals[ii]), 
                              xytext=(8*dpi_scale,-6*dpi_scale), 
                              textcoords='offset pixels',
                            )
        elif names[ii] == 'JET':
            if T_vals[ii] < 5:
                plt.annotate( names[ii], 
                              xy=( T_vals[ii], nTtau_vals[ii]), 
                              xytext=(-14*dpi_scale,6*dpi_scale), 
                              textcoords='offset pixels',
                            )
            elif T_vals[ii] > 25:
                plt.annotate( names[ii], 
                              xy=( T_vals[ii], nTtau_vals[ii]), 
                              xytext=(-10*dpi_scale,-17*dpi_scale), 
                              textcoords='offset pixels',
                            )
        elif names[ii] == 'JT60U':
            if T_vals[ii] > 40:
                plt.annotate( names[ii], 
                              xy=( T_vals[ii], nTtau_vals[ii]), 
                              xytext=(8*dpi_scale,-5*dpi_scale), 
                              textcoords='offset pixels',
                            )
            elif T_vals[ii] < 18:
                # only arrow (to have all arrows starting at same point)
                plt.annotate( '', 
                              xy=( T_vals[ii], nTtau_vals[ii]), 
                              xytext=(11, 1.2e21),
                              arrowprops=dict( arrowstyle="->", shrinkA=0, color='0.2' )
                            )
                # make annotation text (is done once for every case...)
                plt.annotate( names[ii],
                              xy=( T_vals[ii], nTtau_vals[ii]), 
                              xytext=(11, 1.2e21),
                              ha='right', va='bottom'
                            )
        elif names[ii] == 'ASDEX Upgrade':
            # only arrow (to have all arrows starting at same point)
            plt.annotate( '', 
                          xy=( T_vals[ii], nTtau_vals[ii]), 
                          xytext=(17, 5e19),
                          arrowprops=dict( arrowstyle="->", shrinkA=0, color='0.2' )
                        )
            # make annotation text (is done once for every case...)
            plt.annotate( names[ii],
                          xy=( T_vals[ii], nTtau_vals[ii]), 
                          xytext=(17, 5e19),
                          ha='left', va='center'
                        )
        elif names[ii] == 'ASDEX':
            # only arrow (to have all arrows starting at same point)
            plt.annotate( '', 
                          xy=( T_vals[ii], nTtau_vals[ii]), 
                          xytext=(7.5, .9e19),
                          arrowprops=dict( arrowstyle="->", shrinkA=0, color='0.2' )
                        )
            # make annotation text (is done once for every case...)
            plt.annotate( names[ii],
                          xy=( T_vals[ii], nTtau_vals[ii]), 
                          xytext=(7.5, .9e19),
                          ha='left', va='center'
                        )
        elif names[ii] == 'DIII-D':
            # only arrow (to have all arrows starting at same point)
            plt.annotate( '', 
                          xy=( T_vals[ii], nTtau_vals[ii]), 
                          xytext=(8, 5e20),
                          arrowprops=dict( arrowstyle="->", shrinkA=0, color='0.2' )
                        )
            # make annotation text (is done once for every case...)
            plt.annotate( names[ii],
                          xy=( T_vals[ii], nTtau_vals[ii]), 
                          xytext=(8, 5e20),
                          ha='center', va='bottom'
                        )
        elif names[ii] == 'JET (DT)':
            # only arrow (to have all arrows starting at same point)
            plt.annotate( '', 
                          xy=( T_vals[ii], nTtau_vals[ii]), 
                          xytext=(23, 1.4e21),
                          arrowprops=dict( arrowstyle="->", shrinkA=0, color='0.2' )
                        )
            # make annotation text (is done once for every case...)
            plt.annotate( names[ii],
                          xy=( T_vals[ii], nTtau_vals[ii]), 
                          xytext=(23, 1.4e21),
                          ha='center', va='bottom'
                        )

    # format plot
    ax1.set_xlabel( "$T_i$ in keV" )
    ax1.set_ylabel( r"$nT_i\tau_E$ in (keVm$^{-3}$s)" )
    ax1.set_xlim( 0.1, 1e2 )
    ax1.set_ylim( 1e16, 1e22 )
    ax1.set_xscale( 'log' )
    ax1.set_yscale( 'log' )

    ax1.legend( loc='lower right' )

    # write credits into plot, to ensure that people know they can use the plot
    # (somebody once told me, every plot appearing somewhere in the internet
    #  should contain information on how to use it, otherwise it is useless)
    # you probably want to remove it when you make you own plot
    # attribution would still be gratefully acknowledged :)
    # also note that the licence refers only to that specific plot
    # the licence for the code is mentioned above and in the LICENCE file
    credit_str = u'Alf Köhn-Seeman, CC BY-SA 4.0'
    fig1.text( .71, .885, credit_str, fontsize=7 )

    make_plot( fname_plot='triple_product_vs_T.png' )

#;}}}


if __name__ == '__main__':
    main()

