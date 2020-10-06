#!/usr/bin/python
# coding: utf-8
"""
Plot the binding energy per nucleon
"""

__author__      = "Alf Köhn-Seemann"
__email__       = "koehn@igvp.uni-stuttgart.de"
__copyright__   = 'Uni Stuttgart'
__license__     = 'MIT'

# import standard modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from urllib.request import urlopen

# credit string to include at top of plot, to ensure people know they can use the plot
# (someone once told me, every plot appearing somewhere in the internet
#  should contain information on how to use it, otherwise it is useless)
# note that the license refers only to that specific plot
# the license for the code is mentioned in the LICENSE file (and above)
credit_str  = f'{__author__}, CC BY-SA 4.0'

# change default plot formatting
plt.rcParams.update({'font.size': 18})


def read_NIST_data( url='', fname='' ):
#{{{
    """
    Read the atomic weight dataset from NIST from web or file.

    Parameters
    ----------
    url: string
        datasets from NIST can be obtained via queries, here we need the following
        http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some
    fname: string
        alternatively to directly reading the data from NIST, NIST also offers to 
        download an ascii file (which you need to do before running this script)
        
    Returns
    -------
    NIST_dataset: dict
    """

    # dictionary into which the data will be saved
    # note: keys must exactly correspond to identifiers on website / in file
    #       more clever way would be to create the dictionary dynamically
    #       reading (using) the identifiers from the website
    NIST_dataset = { 'Atomic Number':[],
                     'Atomic Symbol':[],
                     'Mass Number':[],
                     'Relative Atomic Mass':[]
                   }

    # if url is provided, read NIST dataset from webpage
    if len(url) > 0:
        # open connection to URL and check if everything is fine
        web_NIST = urlopen(url)
        if (web_NIST.getcode() != 200):
            print( 'ERROR: http-code while trying to read NIST data from web: {0}'.format(web_NIST.getcode()) )
            return -1
        # read website line by line
        for line in web_NIST:
            # decode line into readable string
            decoded_line = line.decode( "utf-8" )

            for label in NIST_dataset:
                if decoded_line.startswith( label ):

                    # get index of "=", starting to look after identifier text at beginning
                    # note: python gets the correct position (text at beginning is just ignored)
                    id_0 = decoded_line.index( "=", len(label) ) + 1

                    # extract value and if present, remove error-margin, given within brackets ()
                    if '(' in decoded_line[ id_0: ]:
                        id_1    = decoded_line.index( "(", len(label) )
                    else:
                        id_1    = len(decoded_line)
                    if label == 'Atomic Symbol':
                        val = decoded_line[id_0:id_1]
                    else:
                        val = float(decoded_line[id_0:id_1])

                    # add value to dictionary
                    NIST_dataset[ label ].append(val)
         

    # if filename is provided during function call, read NIST dataset from file
    if len(fname) > 0:
        # loop through file line by line
        for line in open( fname ):
            if "=" in line:
                tag, value = line.split( '=' )
                tag   = tag.strip()
                value = value.strip()
                # remove parantheses in values which indicate error
                if '(' in value:
                    value = value.split( '(' )[0]
                if tag in NIST_dataset.keys():
                    NIST_dataset[ tag ].append( value )

            # check if list lengths is the same for each key
            if line == '\n':
                lengths = []
                # get the lengths of the arrays for each key
                for key in NIST_dataset:
                    lengths.append( len( NIST_dataset[key] ) )
                # check if all the lengths are equal
                if not all( elem == lengths[0] for elem in lengths ):
                    print( 'WARNING: there is an error in extracting the data for the file' )
                    print( '         and sorting it into a dictionary' )

    return NIST_dataset
#}}}


def get_mass_number( NIST_dataset ):
#{{{
    """
    Extracts the mass number from the NIST dataset.

    Parameters
    ----------
    NIST_dataset: dict

    Returns
    -------
    mass_number: numpy array
    """

    atomic_mass = np.asarray( NIST_dataset['Relative Atomic Mass'], dtype=np.float32 )
    mass_number = np.around( atomic_mass )

    return mass_number
#}}}


def get_binding_energy( NIST_dataset, norm=True):
#{{{
    """

    Parameters
    ----------
    NIST_dataset: dict
    norm: boolean

    Returns
    -------
    binding_energy: numpy array
    """

    # atomic mass unit in units of MeV/c^2
    amu         = 931.494095
    # masses of particle in units of u
    neutron_u   = 1.00866491588
    proton_u    = 1.00727646688
    electron_u  = 0.00054858

    atomic_mass     = np.asarray( NIST_dataset['Relative Atomic Mass'], dtype=np.float32 )
    atomic_number   = np.asarray( NIST_dataset['Atomic Number'], dtype=np.float32 )
    mass_number     = get_mass_number( NIST_dataset )
    n_neutrons      = mass_number - atomic_number
    mass_defect     = ( ( n_neutrons*neutron_u + atomic_number*proton_u + atomic_number*electron_u ) - atomic_mass ) * amu

    if norm:
        bind_energy = mass_defect/mass_number
    else:
        bind_energy = mass_defect

    return bind_energy
#}}}


def main():
#{{{

    # empty string results in plot into window, otherwise into file
    fname_plot  = ''

    # webpage with NIST data
    url = "http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some"
    #fname = 'atomic_weights.dat'

    # load atomic weights dataset from NIST and process data 
    NIST_dataset        = read_NIST_data( url=url )#, #fname=fname )
    mass_number         = get_mass_number( NIST_dataset )
    bind_energy_norm    = get_binding_energy( NIST_dataset, norm=True)

    # start plot
    german_labels = True

    # create empty figure with no axes
    fig = plt.figure( figsize=(10,6) )
    # add a 1x1 grid of axes (note: axes=plots)
    ax = fig.add_subplot(111)

    # plot a vertical line marking the transition from fusion to fission area
    ax.axvline( x=62, color='lightgrey', linewidth=5 )

    # plot the binding energy per nucleon as function of mass number
    ax.plot( mass_number, bind_energy_norm, 'o--', markersize=10 )

    # format the x-axis
    ax.set_xscale( 'log' )
    ax.set_xlim( min(mass_number), max(mass_number) )
    ax.xaxis.set_major_formatter( ScalarFormatter() )
    if german_labels:
        ax.set_xlabel( r'Massenzahl $A$' )
    else:
        ax.set_xlabel( 'atomic mass number A' )

    # format the y-axis
    if german_labels:
        ax.set_ylabel( 'Bindungsenergie pro Nukleon in MeV' )
    else:
        ax.set_ylabel( 'binding energy per nucleon in MeV' )

    # annotate a few important elements
    plt.annotate( '$\mathregular{^1H}$',  
                  xy=(1,0),   xytext=(2.3,.23), size='large',
                  arrowprops=dict(arrowstyle="->" ) )
    plt.annotate( '$\mathregular{^2H}$',  
                  xy=(2,1),    xytext=(1.3,1.3) ,
                  size='large' )
    plt.annotate( '$\mathregular{^3H}$',  
                  xy=(3,2.83), xytext=(2.1,3),
                  size='large' )
    plt.annotate( '$\mathregular{^3He}$',  
                  xy=(3,2.57), xytext=(2.7,1.7) ,
                  size='large' )
    plt.annotate( '$\mathregular{^4He}$', 
            xy=(4,7),   xytext=(3.1,7.3), size='large' )
    plt.annotate( '$\mathregular{^6Li}$', 
            xy=(6,5.3), xytext=(4.5,4.7), size='large' )
    plt.annotate( '$\mathregular{^{9}Be}$', 
            xy=(9,6.5), xytext=(20.,4.), size='large',  
            arrowprops=dict(arrowstyle="->" ) )
    plt.annotate( '$\mathregular{^{11}B}$', 
            xy=(11,6.9), xytext=(30.,6.), size='large',  
            arrowprops=dict(arrowstyle="->" ) )
    plt.annotate( '$\mathregular{^{12}C}$', 
            xy=(12,7.67), xytext=(5.,7.9), size='large',
            arrowprops=dict(arrowstyle="->" ) )
    plt.annotate( '$\mathregular{^{16}O}$', 
            xy=(16,8), xytext=(7.,8.2), size='large',
            arrowprops=dict(arrowstyle="->" ) )
    plt.annotate( '$\mathregular{^{62}Ni}$', 
            xy=(62,8.8), xytext=(35.,7.), size='large',
            arrowprops=dict(arrowstyle="->" ) )
    plt.annotate( '$\mathregular{^{235}U}$', 
                  xy=(235,7.6), xytext=(100.,6.), size='large',
                  arrowprops=dict(arrowstyle="->" ) )

    # indicate area useful for fusion and for fission by arrows
    if german_labels:
        ax.arrow( 13, 1, 20, 0,                             # x0, y0, dx, dy
                  head_width=.5, head_length=7.,
                  edgecolor='black', facecolor='black',
                  label='fusion' )
        ax.annotate( 'Fusion', 
                     xy=(15, 1.2), size='large' )
    else:
        ax.arrow( 13, 1, 20, 0,                             # x0, y0, dx, dy
                  head_width=.5, head_length=5.,
                  edgecolor='black', facecolor='black',
                  label='fusion' )
        ax.annotate( 'fusion', 
                     xy=(15, 1.2), size='large' )
    if german_labels:
        ax.arrow( 260, 1, -165, 0,                          # x0, y0, dx, dy
                  head_width=.5, head_length=15.,
                  edgecolor='black', facecolor='black')
        ax.annotate( 'Spaltung', 
                     xy=(97, 1.2), size='large' )
    else:
        ax.arrow( 240, 1, -140, 0,                          # x0, y0, dx, dy
                  head_width=.5, head_length=15.,
                  edgecolor='black', facecolor='black')
        ax.annotate( 'fission', 
                     xy=(110, 1.2), size='large' )

    # force ticks to point inwards
    ax.tick_params( axis='both', which='both', direction='in', top=True, right=True )

    if len( fname_plot ):
        plt.savefig( fname_plot, bbox_inches='tight', dpi=600 )
        print( '    plot written into {0}'.format( fname_plot ) )
    else:
        plt.show()

#}}}


if __name__ == '__main__':
    main()

