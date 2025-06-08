# coding=utf-8

"""
Plot a field line on a flux surface, where only half of the latter is plotted.

Examples:
    $ python plot_fieldline.py  <-- plot fieldline with default values of q=3.5
    $ python plot_fieldline.py -q 2.7  <-- plot with q=2.7
    $ python plot_fieldline.py -f plot.png  <-- plot into file 'plot.png'
    $ python plot_fieldline.py -n 2  <-- make 2 toroidal circumferences 
                                         (10 is default value)
"""

__author__      = 'Alf Köhn-Seemann'
__email__       = 'koehn@igvp.uni-stuttgart.de'
__copyright__   = 'University of Stuttgart'
__license__     = 'MIT'


# import standard modules
import argparse
import numpy as np
from mayavi import mlab


def toroidal_helix( R0=3., a=1., phi=np.linspace(0., 1.*np.pi, 200), q=2. ):
#{{{

    n = 1./q

    x = (R0 - a*np.cos(n*phi))*np.cos(phi)
    y = (R0 - a*np.cos(n*phi))*np.sin(phi)
    z = a*np.sin(n*phi)

    return x, y, z
#}}}


def make_torus_coordinates( R0, a=1., 
                            theta=np.linspace(0, 2.*np.pi, 200),
                            phi=np.linspace(0, 1.*np.pi, 200)
                          ):
#{{{
    """
    Return x, y, z coordinates of a torus.

    Parameters
    ----------
    R0: float
        major radius of torus
    a: float
        minor radius of torus
    theta: numpy arr (1D) of float
        poloidal angle
    phi: numpy arr (1D) of float
        toroidal angle

    Returns
    -------
    x_torus, y_torus, z_torus: numpy arr (2D) or float
        x-, y-, z-coordinates of a torus as defined by input parameters
    """

    theta, phi  = np.meshgrid(theta, phi)

    x_torus     = (R0 + a*np.cos(theta)) * np.cos(phi)
    y_torus     = (R0 + a*np.cos(theta)) * np.sin(phi)
    z_torus     = a * np.sin(theta)

    return x_torus, y_torus, z_torus
#}}}


def make_plot( R0=3., a=1., fname_plot='', q_fieldline=3.5, fieldline_torTurns=10):
#{{{

    # define coordinates for field line with tor_circumf toroidal circumferences
    phi_fieldline   = np.linspace(0*np.pi, fieldline_torTurns*2.*np.pi, 1000)
    x_fieldline, y_fieldline, z_fieldline = toroidal_helix( R0=R0, a=a, phi=phi_fieldline, q=q_fieldline )

    # define coordinates for flux surfaces (here with circular cross section)
    # note: only plot half a torus, thus phi=0...pi
    theta = np.linspace(0, 2.*np.pi, 200)
    phi   = np.linspace(0, 1.*np.pi, 200)
    # for better visibility of field lines, flux surfaces are made slightly smaller
    x_flufla, y_flufla, z_flufla = make_torus_coordinates( R0, a=a*.99, 
                                                           theta=theta, 
                                                           phi=phi )

    # create the figure and define some basic properties (color, viewangle, size)
    fig1 = mlab.figure( bgcolor=(1,1,1), size=(1000,1000) )
    mlab.view(azimuth=90, elevation=125, figure=fig1)

    # plot fieldline
    mlab.plot3d( x_fieldline, y_fieldline, z_fieldline, 
                 color=(0,0,0), figure=fig1 )

    # plot flux surface
    mlab.mesh( x_flufla, y_flufla, z_flufla, 
               color=(0,.75,1.), figure=fig1 )

    if len(fname_plot) > 0:
        mlab.savefig( fname_plot, size=(4000, 4000) )
    else:
        mlab.show()

#}}}


def main():
#{{{

    # initialize parser for command line options
    parser  = argparse.ArgumentParser()
    # add optional arguments
    parser.add_argument( "-f", "--fname_plot", type=str, default='',
                         help="Filename for plot." )
    parser.add_argument( "-q", "--q_fieldline", type=float, default=3.5, 
                         help="Safety factor of field line." )
    parser.add_argument( "-n", "--nTorTurns", type=float, default=10,
                         help="Toroidal circumferences of field line." )
    # read all arguments from command line
    args    = parser.parse_args()
    fname_plot  = args.fname_plot
    q_fieldline = args.q_fieldline
    tor_circumf = args.nTorTurns

    # major and minor radius
    R0, a = 3., 1.

    make_plot( R0=R0, a=a, fname_plot=fname_plot, q_fieldline=q_fieldline, fieldline_torTurns=tor_circumf )

#}}}


if __name__ == '__main__':
    main()

