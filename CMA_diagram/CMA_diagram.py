# coding=utf-8

"""
Plot a CMA diagram.

In a typical CMA diagram, the horizontal axis, corresponds to X=(omega_pe/omega_0)^2,
where omega_pe is the electron plasma frequency and omega_0 the angular wave frequency,
thus it corresponds to the electron plasma density n_e.
The vertical axis corresponds to Y=(omega_ce/omega_0), where omega_ce is the 
electron cyclotron frequency, thus it corresponds to the absolute value of the 
background magnetic field B_0.

Simply run this script to produce a png plot:
    $ python CMA_diagram.py
"""

__author__      = 'Alf Köhn-Seemann'
__email__       = 'koehn@igvp.uni-stuttgart.de'
__copyright__   = 'University of Stuttgart'
__license__     = 'MIT'

# import standard modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker


plt.rcParams.update({'font.size':12})
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


def oplot_Ocut( ax, y_range=[], linestyle='solid', linewidth=3, color='black' ):
    #;{{{
    """
    Overplot the O-mode cut-off in a CMA diagram.

    In a typical CMA diagram, the horizontal axis corresponds to X=(omega_pe/omega_0)^2,
    where omega_pe is the electron plasma frequency and omega_0 the angular wave frequency;
    the vertical axis corresponds to Y=(omega_ce/omega_0), where omega_ce is the 
    electron cyclotron frequency.
    The O-mode cut-off is thus located at X=1.

    Parameters
    ----------
    ax: Axes object
    y_range: list or np.array
        2-element list or array specifying the start and end point for the line 
        indicating the O-mode cut-off. If not provided, the range of the y-axis
        is used.
    linestyle: str
    linewidth: int
    color: str

    Returns
    -------
    """
    if len(y_range) == 0:
        y0, y1 = ax.get_ylim()
    else:
        y0  = y_range[0]
        y1  = y_range[1]

    # plot the cut-off position
    ax.plot( [1,1], [y0,y1],
             marker='None',
             linestyle=linestyle, linewidth=linewidth, 
             color=color, 
           )
    # write text to the cut-off (annotate it)
    txt_y   = y1 - .3*(y1-y0)   #1.4
    ax.annotate( 'O cut-off', xy=(1.,txt_y), xytext=(1.,txt_y), rotation=90,
                 horizontalalignment='right', verticalalignment='bottom',
               )
    #;}}}


def oplot_XRcut( ax, linestyle='solid', linewidth=3, color='black' ):
    #;{{{
    """
    Overplot the X-mode R cut-off in a CMA diagram.

    In a typical CMA diagram, the horizontal axis corresponds to X=(omega_pe/omega_0)^2,
    where omega_pe is the electron plasma frequency and omega_0 the angular wave frequency;
    the vertical axis corresponds to Y=(omega_ce/omega_0), where omega_ce is the 
    electron cyclotron frequency.
    The right-hand cut-off reads
        w_Rcut = 0.5*( sqrt(w_ce^2+4*w_pe^2) + w_ce )
        <=>
        w_Rcut - 0.5*w_ce = 0.5 * sqrt(w_ce^2+4*w_pe^2)
    Normalizing to w_0 yields
        2 - Y = sqrt(Y^2 + 4*X)
        <=>
        4 - 4Y + Y^2 = Y^2 + 4*X
        <=>
        4*(1-Y) = 4*X
        <=>
        1-Y = X
        <=>
        Y = 1 - X

    Parameters
    ----------
    ax: Axes object
    linestyle: str
    linewidth: int
    color: str

    Returns
    ------
    """

    x0, x1  = ax.get_ylim()
    x_range = np.array( [x0, x1] )

    arr_x   = np.linspace( np.min(x_range), np.max(x_range), 200 )

    # calculate right-hand cut-off
    arr_y   = 1. - arr_x 
    ax.plot( arr_x, arr_y, 
             linestyle=linestyle, marker='None', linewidth=linewidth, 
             color=color, 
           )

    ax.annotate( 'XR cut-off', xy=(.2,.2), xytext=(.2,.2), rotation=-47,
                  horizontalalignment='left', verticalalignment='bottom',
               )
    #;}}}


def oplot_XLcut( ax, x_range=[], linestyle='solid', linewidth=3, color='black' ):
    #;{{{
    """
    Overplot the X-mode L cut-off in a CMA diagram.

    In a typical CMA diagram, the horizontal axis corresponds to X=(omega_pe/omega_0)^2,
    where omega_pe is the electron plasma frequency and omega_0 the angular wave frequency;
    the vertical axis corresponds to Y=(omega_ce/omega_0), where omega_ce is the 
    electron cyclotron frequency.
    The left-hand cut-off reads
        w_Lcut = 0.5*( sqrt(w_ce^2+4*w_pe^2) - w_ce )
        <=>
        w_Rcut + 0.5*w_ce = 0.5 * sqrt(w_ce^2+4*w_pe^2)
    Normalizing to w_0 yields
        2 + Y = sqrt(Y^2 + 4*X)
        <=>
        4 + 4Y + Y^2 = Y^2 + 4*X
        <=>
        4*(1+Y) = 4*X
        <=>
        1+Y = X
        <=>
        Y = X - 1

    Parameters
    ----------
    ax: Axes object
    x_range: list or np.array
        2-element list or array specifying the start and end point for the line 
        indicating the X-mode L cut-off. If not provided, the range of the x-axis
        is used.
    linestyle: str
    linewidth: int
    color: str

    Returns
    ------
    """

    if len(x_range) == 0:
        x0, x1  = ax.get_xlim()
        x_range = np.array( [x0, x1] )

    arr_x = np.linspace( 1., np.max(x_range), 200 )
    arr_y = -1. + arr_x
    ax.plot( arr_x, arr_y,
             marker='None',
             linestyle=linestyle, linewidth=linewidth, 
             color=color, 
            )
    ax.annotate( 'XL cut-off', xy=(1.2,.33), xytext=(1.2,.33), rotation=47,
                 horizontalalignment='left', verticalalignment='bottom',
               )
    #;}}}


def oplot_Xres( ax, theta=np.array([90.]), annotation_x=np.array([.5]), 
                linestyle='dashed', linewidth=3, color='black' 
              ):
    #;{{{
    """
    Overplot the X-mode upper-hybrid resonance in a CMA diagram.

    In a typical CMA diagram, the horizontal axis corresponds to X=(omega_pe/omega_0)^2,
    where omega_pe is the electron plasma frequency and omega_0 the angular wave frequency;
    the vertical axis corresponds to Y=(omega_ce/omega_0), where omega_ce is the 
    electron cyclotron frequency.
    
    Parameters
    ----------
    ax: Axes object
    theta: np.array
        Angle of microwave propagation with respect to the background magnetic field B_0,
        units are degrees.
    annotation_x: np.array
        X-coordinates for the labels (which are the theta-values) to be written into the plot.
    linestyle: str
    linewidth: int
    color: str

    Returns
    ------
    """

    for ii in range( len(theta) ):

        arr_x = np.linspace( 0., 1., 200 )
        arr_y = +1.*np.sqrt( (1.-arr_x)/(1.-arr_x*np.cos(theta[ii]/180.*np.pi))  )

        # remove NaN-values in arr_y
        arr_x = arr_x[ np.isfinite(arr_y) ]
        arr_y = arr_y[ np.isfinite(arr_y) ]

        # string for legend (only for one theta value)
        if theta[ii] == theta[0]:
            label_str = 'X-resonance'
        else:
            label_str = None
        ax.plot( arr_x, arr_y,
                 marker='None',
                 linestyle=linestyle, linewidth=linewidth,
                 color=color, 
                 label=label_str
                )

        # annotate with theta-value
        annot_x_id = np.argmin( np.abs(arr_x-annotation_x[ii]) )
        annot_x = arr_x[annot_x_id]
        annot_y = arr_y[annot_x_id]
        annot_txt = r'${0:2.0f}\degree$'.format( theta[ii] )
        ax.annotate( annot_txt, xy=(annot_x,annot_y), xytext=(annot_x,annot_y-.05), 
                      horizontalalignment='left', verticalalignment='top',
                    )
    #;}}}


def oplot_Ores( ax, theta=np.array([30.]), annotation_x=np.array([1.5]),
                x_range=[], 
                linestyle='dotted', linewidth=3, color='black' 
              ):
    #;{{{
    """
    Overplot the O-mode O resonance in a CMA diagram.

    In a typical CMA diagram, the horizontal axis corresponds to X=(omega_pe/omega_0)^2,
    where omega_pe is the electron plasma frequency and omega_0 the angular wave frequency;
    the vertical axis corresponds to Y=(omega_ce/omega_0), where omega_ce is the 
    electron cyclotron frequency.
    
    Parameters
    ----------
    ax: Axes object
    theta: np.array
        Angle of microwave propagation with respect to the background magnetic field B_0,
        units are degrees.
    annotation_x: np.array
        X-coordinates for the labels (which are the theta-values) to be written into the plot.
    x_range: list or np.array
        2-element list or array, where currently only the second element is used to specify 
        the end point for the line indicating the O-resonance. If not provided, the range of 
        the x-axis is used.
    linestyle: str
    linewidth: int
    color: str

    Returns
    """

    if len(x_range) == 0:
        x0, x1  = ax.get_xlim()
        x_range = np.array( [x0, x1] )

    arr_x = np.linspace( (1.+1e-6), np.max(x_range), 200)

    for ii in range( len(theta) ):
        arr_y = +1.*np.sqrt( (1.-arr_x)/(1.-arr_x*np.cos(theta[ii]/180.*np.pi))  )

        # string for legend (only for one theta value)
        if theta[ii] == theta[0]:
            label_str = 'O-resonance'
        else:
            label_str = None
        ax.plot( arr_x, arr_y,
                 marker='None',
                 linestyle=linestyle, linewidth=linewidth, 
                 color=color, 
                 label=label_str
               )

        # annotate with theta-value
        annot_x_id = np.argmin( np.abs(arr_x-annotation_x[ii]) )
        annot_x = arr_x[annot_x_id]
        annot_y = arr_y[annot_x_id]
        annot_txt = r'${0:2.0f}\degree$'.format( theta[ii] )
        ax.annotate( annot_txt, xy=(annot_x,annot_y), xytext=(annot_x,annot_y+.01), 
                      horizontalalignment='left', verticalalignment='bottom',
                    )

    #;}}}


def oplot_ECR( ax, x_range=[], linestyle='solid', linewidth=3, color='black' ):
    #;{{{
    """
    Overplot the electron cyclotron resonance in a CMA diagram.

    In a typical CMA diagram, the horizontal axis corresponds to X=(omega_pe/omega_0)^2,
    where omega_pe is the electron plasma frequency and omega_0 the angular wave frequency;
    the vertical axis corresponds to Y=(omega_ce/omega_0), where omega_ce is the 

    Parameters
    ----------
    ax: Axes object
    x_range: list or np.array
        2-element list or array specifying the start and end point for the line 
        indicating the electron cyclotron resonance. If not provided, the range of the x-axis
        is used.
    linestyle: str
    linewidth: int
    color: str

    Returns
    ------
    """

    if len(x_range) == 0:
        x0, x1  = ax.get_xlim()
        x_range = np.array( [x0, x1] )

    arr_x = x_range
    arr_y = np.array( [1.,1.] )
    ax.plot( arr_x, arr_y,
             marker='None', 
             linestyle=linestyle, linewidth=linewidth,
             color=color, 
            )
    label_str_Rres = 'ECR'
    ax.annotate( label_str_Rres, xy=(.1,1), xytext=(.1,1.02),
                 horizontalalignment='left', verticalalignment='bottom',     
               )
    #;}}}
                

def main():
    #;{{{

    # plot configuration
    fname_plot  = 'CMA_diagram.png'
    # linewidth for O- and X-mode
    lw_O        = 3
    lw_X        = 3
    # color for O- and X-mode
    color_O     = 'black'
    color_X     = 'black' #'red'
    # linestyle for cut-offs and resonances
    ls_Ocut     = 'solid'
    ls_Xcut     = 'solid'
    ls_Xres     = 'dashed'
    ls_Ores     = 'dotted'

    # (width, height) in inches
    fig = plt.figure( figsize=(8,6) )

    ax1 = fig.add_subplot( 1,1,1 )

    x_range = np.array( [0,3] )
    y_range = np.array( [0,2] )

    # oplot O cut-off
    oplot_Ocut( ax1, y_range=y_range, linestyle=ls_Ocut, linewidth=lw_O, color=color_O )

    # oplot XR cut-off
    oplot_XRcut( ax1, linestyle=ls_Xcut, linewidth=lw_X, color=color_X )

    # oplot XL cut-off
    oplot_XLcut( ax1, x_range=x_range, linestyle=ls_Xcut, linewidth=lw_X, color=color_X )

    # oplot resonances for different thetas
    thetas          = np.array( [90., 30., 10.] )
    annotations_x   = np.array( [.5, .65, .8] )
    # X resonance
    oplot_Xres( ax1, theta=thetas, annotation_x=annotations_x, 
                linestyle=ls_Xres, linewidth=lw_X, color=color_X 
              )
    # O resonance
    thetas          = np.array( [30., 10.] )
    annotations_x   = np.array( [1.5, 1.2] )
    oplot_Ores( ax1, theta=thetas, annotation_x=annotations_x,
                x_range=x_range,
                linestyle=ls_Ores, linewidth=lw_O, color=color_O 
              )
    # R-resonance (resonance for theta=0)
    oplot_ECR( ax1, x_range=x_range, linestyle='solid', linewidth=lw_X, color=color_X )

    ax1.legend( loc='lower right' )

    # set axes labels
    ax1.set_xlabel( r'$X = (\omega_{pe}/\omega_0)^2 \propto n_e$' )
    ax1.set_ylabel( r'$Y = \omega_{ce}/\omega_0 \propto |B_0|$' )

    # set plot ranges
    ax1.set_xlim( x_range )
    ax1.set_ylim( y_range )

    ax1.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter("{x:.1f}"))
    ax1.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter("{x:.1f}"))

    # write credits into plot, to ensure that people know they can use the plot
    # (somebody once told me, every plot appearing somewhere in the internet
    #  should contain information on how to use it, otherwise it is useless)
    # you probably want to remove it when you make you own plot
    # attribution would still be gratefully acknowledged :)
    # also note that the licence refers only to that specific plot
    # the licence for the code is mentioned above and in the LICENCE file
    credit_str = u'Alf Köhn-Seeman, CC BY-SA 4.0'
    fig.text( .71, .885, credit_str, fontsize=7 )

    make_plot( fname_plot=fname_plot )
    #;}}}


if __name__ == '__main__':
    main()

