from datetime import date
from math import atan2, degrees

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

from matplotlib import container
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerErrorbar

from matplotlib.path import Path
from matplotlib.patches import PathPatch

import colorsys
import matplotlib.colors as mcolors
import matplotlib.cm as cm


# ---------------------------------------------------
# Formatting:
# ---------------------------------------------------
_small_size = 10
_medium_size = 12
_bigger_size = 14
_large_size = 16

plt.rc('font', size=_medium_size)         # controls default text sizes
plt.rc('axes', titlesize=_bigger_size)    # fontsize of the axes title
plt.rc('axes', labelsize=_bigger_size)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=_small_size)    # fontsize of the tick labels
plt.rc('ytick', labelsize=_small_size)    # fontsize of the tick labels
plt.rc('legend', fontsize=_medium_size)   # legend fontsize
plt.rc('figure', titlesize=_large_size)   # fontsize of the figure title

# ---------------------------------------------------
# Styles:
# ---------------------------------------------------
# Line plot style
style_solid = {'ls':'-', 'lw':2}
style_dashed = {'ls':'--', 'lw':1.5}

# Scatter plot style
style_scatter = {'s':2.0, 'alpha':.9}

# Errorbar plot styles
style_yerr = {'xerr':0, 'markersize':2.5, 'fmt':'s',
              'elinewidth':2.3, 'capsize':0, 'zorder':1}
style_yerr_ps = {'xerr':0, 'markersize':4, 'fmt':'o',
              'elinewidth':2.3, 'capsize':0, 'zorder':.5}
# MIT open data style (P. Komiske)
modstyle = {'lw':2, 'capsize':2, 'capthick':1.5, 'markersize':2,
            'linestyle':'None', 'zorder':1}
modstyle_ps = {'lw':2, 'capsize':2, 'capthick':1.5, 'markersize':5,
               'fmt':'o', 'linestyle':'None', 'zorder':.5}

# ####################################
# Plot creation
# ####################################
# ---------------------------------------------------
# Basic figure type:
# ---------------------------------------------------

def aestheticfig(xlabel='x', ylabel=r'Probability Density',
                 title=None, showdate=True,
                 xlim=(0, 1), ylim=(0, 1), ylim_ratio=(0.5, 2.),
                 ratio_plot=True, ylabel_ratio='Ratio',
                 labeltext=None):
    """Creates a figure and associated axes. Can be used to
    produce a figure with a subplot which is, for example,
    associated with a ratio.

    Parameters
    ----------
    xlabel : str
        xlabel of the plot.
    ylabel : str
        ylabel of the plot.
    title : str
        title of the plot.
    showdate : bool
        If True, adds a date to the upper right of the plot.
    xlim : tuple
        The x limits of the plot.
    ylim : tuple
        The y limits of the plot.
    ylim_ratio : tuple
        The y limits of the ratio subplot.
    ratio_plot : bool
        Determines whether there is an additional subplot
        for ratio plotting.
    ylabel_ratio : str
        ylabel of the ratio subplot, if it exists.

    Returns
    -------
    Figure, axes.Axes
        The figure and axes/subplots specified by the
        above parameters.
    """
    # aesthetic options
    # fig_width = 5.
    # golden_mean = (np.sqrt(5)-1.0)/2.0
    # fig_height = fig_width/golden_mean
    fig_width = 6.4
    fig_height = 4.8
    figsize = (fig_width, fig_height)

    gridspec_kw = {'height_ratios': (3.5, 1) if ratio_plot else (1,),
                   'hspace': 0.0}

    # get subplots
    nsubplots = 2 if ratio_plot else 1
    fig, axes = plt.subplots(nsubplots, gridspec_kw=gridspec_kw,
                             figsize=figsize)
    if nsubplots == 1:
        axes = [axes]

    # axes limits
    for ax in axes:
        ax.set_xlim(*xlim)
    axes[0].set_ylim(*ylim)
    if ratio_plot:
        axes[1].set_ylim(*ylim_ratio)
        axes[1].set_yscale('log')

    # axes labels
    axes[-1].set_xlabel(xlabel)
    axes[0].set_ylabel(ylabel, labelpad=5)
    if ratio_plot:
        axes[1].set_ylabel(ylabel_ratio, labelpad=-10)

    # tick settings
    for ax in axes:
        ax.minorticks_on()
        ax.tick_params(top=True, right=True, bottom=True,
                       left=True, direction='in', which='both')

    if ratio_plot:
        axes[0].tick_params(labelbottom=False)
        axes[1].tick_params(axis='y')

    # Extra plot information
    pad = .01

    if showdate:
        # Including date
        axes[0].text(
            x=1,
            y=1.005+pad,
            s=date.today().strftime("%m/%d/%y"),
            transform=axes[0].transAxes,
            ha="right",
            va="bottom",
            fontsize=_medium_size * 0.95,
            fontweight="normal"
        )

    if labeltext is not None:
        # Extra primary label
        axes[0].text(
            x=-0.1,
            y=1.005+pad,
            s=labeltext,
            transform=axes[0].transAxes,
            ha="left",
            va="bottom",
            fontsize=_medium_size * 1.5,
            fontweight="bold",
            fontname="DIN Condensed"
        )

    if title is not None:
        # Main title
        axes[0].text(
            x=.12,
            y=1.005+pad,
            s=title,
            transform=axes[0].transAxes,
            ha="left",
            va="bottom",
            fontsize=_medium_size * 1.5,
            fontstyle="italic",
            fontname="Arial"
        )

    plt.tight_layout()

    return fig, axes


def aesthetic_N_by_M(xlabel='x', ylabel=r'Probability Density',
                 nrows=2, ncols=2,
                 title=None, showdate=False,
                 xlim=(0, 1), ylim=(0, 1), ylim_ratio=(0.5, 2.),
                 ratio_plot=True, ylabel_ratio='Ratio',
                 labeltext=None,
                 colnames=None,
                 rownames=None):
    """Creates a figure and associated axes. Can be used to
    produce a figure with a subplot which is, for example,
    associated with a ratio.

    Parameters
    ----------
    xlabel : str
        xlabel of the plot.
    ylabel : str
        ylabel of the plot.
    nrows : int
        number of rows in the plot.
    ncols : int
        number of columns in the plot.
    title : str
        title of the plot.
    showdate : bool
        If True, adds a date to the upper right of the plot.
    xlim : tuple
        The x limits of the plot.
    ylim : tuple
        The y limits of the plot.
    ylim_ratio : tuple
        The y limits of the ratio subplot.
    ratio_plot : bool
        Determines whether there is an additional subplot
        for ratio plotting.
    ylabel_ratio : str
        ylabel of the ratio subplot, if it exists.
    colnames : list of str
        labels of the columns for the grid of plots.
    rownames : list of str
        labels of the rows for the grid of plots.

    Returns
    -------
    Figure, axes.Axes
        The figure and axes/subplots specified by the
        above parameters.
    """
    # aesthetic options
    # fig_width = 5.
    # golden_mean = (np.sqrt(5)-1.0)/2.0
    # fig_height = fig_width/golden_mean
    fig_width = 6.4 * 1.5
    fig_height = 4.8 * 1.5
    figsize = (fig_width, fig_height)

    nsubplots = 2 if ratio_plot else 1

    gridspec_kw = {'height_ratios': (3.5, 1) if ratio_plot else (1,),
                   'hspace': 0.0}

    # Get figure
    fig = plt.figure(figsize=figsize)
    axes = []

    # Set up grid of subplots
    gs_full = GridSpec(nrows, ncols, wspace=0.3, hspace=0.3)

    for irow in range(nrows):
        for jcol in range(ncols):
            if ratio_plot:
               gs_ratio=GridSpecFromSubplotSpec(2, 1,
                                 hspace=0.0,
                                 height_ratios=(3.5, 1),
                                 subplot_spec=gs_full[irow,jcol])
               ax_main = fig.add_subplot(gs_ratio[0])
               ax_ratio = fig.add_subplot(gs_ratio[1])
               axes.append(ax_main)
               axes.append(ax_ratio)
            else:
                axes.append(fig.add_subplot(gs_full[irow,jcol]))

    # reesetting indices in case we have row and column labels
    irow, jcol = 0, 0

    # Setting up Axes 
    for i, ax in enumerate(axes):
        # ---------------------------------
        # Universal settings
        # ---------------------------------
        # x- and y-limits 
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)

        # Ticks
        ax.minorticks_on()
        ax.tick_params(top=True, right=True, bottom=True,
                       left=True, direction='in', which='both')
        
        # ---------------------------------
        # Local settings
        # ---------------------------------
        # Location properties of different axes
        is_ratio = (ratio_plot and i%2 == 1)
        is_leftmost = (i%nrows == 0) if not ratio_plot\
            else (i%(2*nrows) in [0,1])
        is_bottom = (i >= len(axes) - nrows) if not ratio_plot\
            else (i >= len(axes) - 2.*nrows and is_ratio)
        is_top = (i < nrows) if not ratio_plot\
            else (i < 2*nrows and not is_ratio)

        # Ratio plots have different y-axis settings 
        if is_ratio:
            ax.set_ylim(ylim_ratio)
            ax.set_yscale('log')
            ax.tick_params(axis='y')
        elif ratio_plot and not is_ratio:
            ax.tick_params(labelbottom=False)

        # Lowest level plots get x-axis labels
        if is_bottom:
            ax.set_xlabel(xlabel)

        # Leftmost plots get y-axis labels
        if is_leftmost:
            ax.set_ylabel(ylabel_ratio if is_ratio else ylabel,
                          labelpad=0 if is_ratio else 5,
                          fontsize=_bigger_size if not is_ratio\
                            else _small_size)
            if not is_ratio and rownames is not None:
                ax.annotate(rownames[irow], xy=(0, 0.5),
                        xytext=(-ax.yaxis.labelpad - 10, 0),
                        xycoords=ax.yaxis.label, textcoords='offset points',
                        size='large', ha='right', va='center',
                        rotation=90)
                irow+=1


        if is_top and colnames is not None:
            ax.annotate(colnames[jcol], xy=(0.5, 1), xytext=(0, 5),
                        xycoords='axes fraction', textcoords='offset points',
                        size='large', ha='center', va='baseline')
            jcol+=1

    # Extra plot information
    pad = .01

    if showdate:
        # Including date
        axes[0].text(
            x=1,
            y=1.005+pad,
            s=date.today().strftime("%m/%d/%y"),
            transform=axes[0].transAxes,
            ha="right",
            va="bottom",
            fontsize=_medium_size * 0.95,
            fontweight="normal"
        )

    if labeltext is not None:
        # Extra primary label
        axes[0].text(
            x=-0.2,
            y=1.005+pad,
            s=labeltext,
            transform=axes[0].transAxes,
            ha="left",
            va="bottom",
            fontsize=_medium_size * 1.5,
            fontweight="bold",
            fontname="DIN Condensed"
        )

    if title is not None:
        # Main title
        axes[0].text(
            x=1.1,
            y=1.15+pad,
            s=title,
            transform=axes[0].transAxes,
            ha="center",
            va="bottom",
            fontsize=_medium_size * 1.5,
            fontstyle="italic",
            fontname="Arial"
        )

    plt.tight_layout()

    if rownames is not None:
        fig.subplots_adjust(left=0.15)
    # if colnames is not None:
    #     fig.subplots_adjust(top=0.95)

    return fig, axes


# ---------------------------------
# Histograms
# ---------------------------------

def get_hist_edges_centers(obs_vals, weights, nbins, binspace,
                           lbin=None, rbin=None):
    """
    Takes in a set of values and an associated binning space.
    Returns (histogram heights, bin edges, bin centers)
    associated with those values.

    Parameters
    ----------
        obs_vals :  The set of values to histogram.
        nbins :     The number of desired bins.
        binspace :  The desired bin space ('lin'/'linear' or 'log').
        lbin, rbin: The left and right edges of the left and rightmost bins.
                    If None, determined by the minimum or maximum
                    values of the given observables, respectively.

    Returns
    -------
        (hist, bin_edges, bin_centers) :
            Three lists describing the properties of the histogram.

    """
    lin_labels, log_labels = ['linear', 'lin'], ['log']

    o_vals = obs_vals if binspace in lin_labels\
             else np.log(obs_vals) if binspace in log_labels\
             else None
    
    good_inds = [np.isfinite(w) and np.isfinite(o)
                 for w, o in zip(weights, o_vals)]
    weights = weights[good_inds]
    o_vals = o_vals[good_inds]

    lbin = lbin if binspace in lin_labels or lbin is None\
             else np.log(lbin)
    rbin = rbin if binspace in lin_labels or rbin is None\
             else np.log(rbin)

    # Making and storing weighted and normalized histogram
    bin_edges = np.linspace(np.nanmin(o_vals) if lbin is None else lbin,
                            np.nanmax(o_vals) if rbin is None else rbin,
                            nbins+1)

    hist, _ = np.histogram(o_vals, bin_edges, weights=weights,
                            density=True)
    if hist[0] == np.nan:
        raise AssertionError

    if binspace in ['linear', 'lin']:
        bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.

    elif binspace in ['log']:
        # Getting logarithmic hists for appropriate observables
        bin_centers = np.exp((bin_edges[:-1] +
                              bin_edges[1:]) / 2.)
        bin_edges = np.exp(bin_edges)

    return (hist, bin_edges, bin_centers)


def text_to_list(filename, use_cols=None, use_rows=None,
                 *kwarg_keys):
    # - - - - - - - - - - - - - - - - -
    # Reading file
    # - - - - - - - - - - - - - - - - -
    data = []

    # Allowing us to store and return extra data from the file,
    if kwarg_keys != ():  # if we are given extra info to look for
        kwargs = {key: None for key in kwarg_keys}

    with open(filename, "r") as file:
        for irow, line in enumerate(file):
            # Read and format the line
            line = line.rstrip("\n")
            info = line.split(" ")

            # Additional information:
            if info[0] == "#" or '' in info:
                # Allowing commented or empty lines
                continue
            if info[0] in kwarg_keys:
                # Assuming we have a row of the form
                # ```[key] = [value]```
                # in the given text file
                kwargs[info[0]] = info[2]
                continue

            # Selecting rows based on given criterion
            if use_rows is not None:
                # If we are explicitly given a list of rows
                if isinstance(use_rows, list) or\
                 isinstance(use_rows, np.ndarray):
                    if irow not in use_rows:
                        continue
                # If we are given a function on the row info
                elif not use_rows(info):
                    continue

            # Default: use all data 
            if use_cols is None:
                data.append(np.array(info))
            # Otherwise, if use_cols is a value or a [value]
            elif (not isinstance(use_cols, list) and\
                  not isinstance(use_cols, np.ndarray)):
                data.append(info[use_cols])
            elif len(use_cols) == 1:
                data.append(info[use_cols[0]])
            # Otherwise, if use_cols is a list
            else:
                data.append([val for i, val in enumerate(info)
                             if i in use_cols])

    if kwarg_keys != ():
        return np.array(data), kwargs
    return np.array(data)


def text_to_hist(filename, use_cols=0, use_rows=None,
                 col_func=None,
                 fig=None, ax=None,
                 nbins=100, binspace='lin',
                 weight_fn=lambda x: np.ones(len(x)),
                 plot_style='plot', save=None,
                 color=None,
                 **kwargs):
    # Default arguments
    if color is None:
        color='cornflowerblue'

    # Getting data from text file
    data = text_to_list(filename, use_cols=use_cols,
                        use_rows=use_rows)

    # If we receive a list of columns to use
    if isinstance(use_cols, list) or isinstance(use_cols, np.ndarray):
        if len(use_cols) != 1:
            # Turning multiple column values into a single number to hist
            assert col_func is not None, "Missing a required function"\
              +" on the columns of the text file for making a histogram."
            data = np.array([col_func(*datum) for datum in data])

    data = data.astype(np.float64)

    # - - - - - - - - - - - - - - - - -
    # Making histogram
    # - - - - - - - - - - - - - - - - -
    # Set up figure, axis limits
    if fig is None and ax is None: 
        fig, ax = aestheticfig(ratio_plot=False, showdate=False, **kwargs)

    ys, x_edges, xs = get_hist_edges_centers(data, weights=weight_fn(data),
                                             nbins=nbins,
                                             binspace=binspace)
    if plot_style == 'plot':
        ax[0].plot(xs, ys, color=color, **style_solid)
    elif plot_style == 'errorbar':
        ax[0].errorbar(xs, ys, yerr=np.sqrt(ys/len(ys)),
                   xerr=(xs - x_edges[:-1], x_edges[1:] - xs),
                   color=color, **modstyle)
    elif plot_style == 'scatter':
        ax[0].scatter(xs, ys, color=color, **style_scatter)

    ax[0].autoscale()
    if binspace == 'log':
        ax[0].set_xscale('log')
        ax[0].set_ylabel('Log-normalized Probability Density')
    else:
        plt.xlim((0, 50))
        plt.ylim(bottom=0)

    if save is not None:
        fig.savefig(save, bbox_inches='tight')
    return


# ---------------------------------------------------
# Putting text on figures:
# ---------------------------------------------------

def set_figtext(fig, text, loc, rightjustify=False, color='black'):
    """Puts text of a nice style on figures."""
    if rightjustify:
        ha = 'right'
    else:
        ha = 'left'

    t = fig.text(*loc, text, linespacing=1.5, ha=ha, color=color)
    t.set_bbox(dict(facecolor='white', alpha=0.9,
                    edgecolor='lightgrey',
                    boxstyle="round,pad=0.35"))

# function to add a stamp to figures
def stamp(left_x, top_y, ax,
          delta_y=0.075,
          textops_update=None,
          **kwargs):

    # text options
    textops = {'horizontalalignment': 'left',
               'verticalalignment': 'center',
               'fontsize': 8.5,
               'transform': ax.transAxes}
    if isinstance(textops_update, dict):
        textops.update(textops_update)

    # add text line by line
    for i in range(len(kwargs)):
        y = top_y - i*delta_y
        t = kwargs.get('line_' + str(i))
        if t is not None:
            ax.text(left_x, y, t, **textops)


# ####################################
# Colors 
# ####################################

def adjust_lightness(color, amount=0.5):
    """
    Adjusts the lightness of the given color by multiplying
    (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> adjust_lightness('g', 0.3)
    >> adjust_lightness('#F034A3', 0.6)
    >> adjust_lightness((.3,.55,.1), 0.5)

    From https://stackoverflow.com/a/49601444
    """
    try:
        c = mcolors.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mcolors.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


def get_colors_colorbar(vals):
    """
    Sets up a colorbar for a set of values to be plotted;
    See, e.g., https://stackoverflow.com/a/30781043.

    Parameters
    ----------
        vals : Values to be plotted

    Returns
    -------
        Colors for each value and a 

    """
    # See, e.g., https://stackoverflow.com/a/30781043
    # setup the normalization and the colormap
    normalize = mcolors.Normalize(vmin=np.min(vals), vmax=np.max(vals))
    colormap = cm.jet

    # List of colors
    colors = [colormap(normalize(rad)) for rad in vals]

    # Associated colorbar (usage: plt.colorbar(scalarmappaple))
    scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
    scalarmappaple.set_array(vals)

    return colors, scalarmappaple


# ####################################
# Error Bands
# ####################################
# From https://matplotlib.org/stable/gallery/lines_bars_and_markers/curve_error_band.html
def draw_error_band(ax, x, y, err, normal=False, **kwargs):
    if normal:
        # Calculate normals via centered finite differences (except the first point
        # which uses a forward difference and the last point which uses a backward
        # difference).
        dx = np.concatenate([[x[1] - x[0]], x[2:] - x[:-2], [x[-1] - x[-2]]])
        dy = np.concatenate([[y[1] - y[0]], y[2:] - y[:-2], [y[-1] - y[-2]]])
        l = np.hypot(dx, dy)
        nx = dy / l
        ny = -dx / l
        
        # end points of errors
        xp = x + nx * err
        yp = y + ny * err
        xn = x - nx * err
        yn = y - ny * err
    else:
        xn, xp = x, x
        yn, yp = y - err, y + err

    vertices = np.block([[xp, xn[::-1]],
                         [yp, yn[::-1]]]).T
    codes = np.full(len(vertices), Path.LINETO)
    codes[0] = codes[len(xp)] = Path.MOVETO
    path = Path(vertices, codes)
    ax.add_patch(PathPatch(path, **kwargs))


# ####################################
# Labelling lines
# ####################################
# From https://stackoverflow.com/a/39402483

def labelLine(line, x, label=None, align=False, **kwargs):
    """Labels a line with the corresponding
    line2d label data.
    """
    ax = line.axes
    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if (x < xdata[0]) or (x > xdata[-1]):
        print('x label location is outside data range!')
        return

    #Find corresponding y co-ordinate and angle of the line
    ip = 1
    for i, xdata_i in enumerate(xdata):
        if x < xdata_i:
            ip = i
            break

    y = (ydata[ip-1] +
         (ydata[ip]-ydata[ip-1])
         *(x-xdata[ip-1])
         /(xdata[ip]-xdata[ip-1]))

    if not label:
        label = line.get_label()

    if align:
        #Compute the slope
        dx = xdata[ip] - xdata[ip-1]
        dy = ydata[ip] - ydata[ip-1]
        ang = degrees(atan2(dy, dx))

        #Transform to screen co-ordinates
        pt = np.array([x, y]).reshape((1, 2))
        trans_angle = ax.transData.transform_angles(
            np.array((ang,)), pt)[0]

    else:
        trans_angle = 0

    #Set a bunch of keyword arguments
    if 'color' not in kwargs:
        kwargs['color'] = line.get_color()

    if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
        kwargs['ha'] = 'center'

    if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
        kwargs['va'] = 'center'

    if 'backgroundcolor' not in kwargs:
        kwargs['backgroundcolor'] = ax.get_facecolor()

    if 'clip_on' not in kwargs:
        kwargs['clip_on'] = True

    if 'zorder' not in kwargs:
        kwargs['zorder'] = 2.5

    t = ax.text(x, y, label, rotation=trans_angle, **kwargs)
    t.set_bbox(dict(facecolor='white', alpha=0.9,
                    edgecolor='lightgrey',
                    boxstyle="round,pad=0.15"))


def labelLines(lines, align=False, xvals=None,
               spacing=None, **kwargs):
    """Labels a set of lines with the
    corresponding line2D label data.
    """
    ax = lines[0].axes
    labLines = []
    labels = []

    #Take only the lines which have labels other than the default ones
    for line in lines:
        label = line.get_label()
        if "_line" not in label:
            labLines.append(line)
            labels.append(label)

    if xvals is None:
        xmin, xmax = ax.get_xlim()
        if spacing is None or spacing == 'lin':
            xvals = np.linspace(xmin, xmax, len(labLines)+2)[1:-1]
        elif spacing == 'log':
            xvals = np.logspace(np.log10(xmin), np.log10(xmax),
                                len(labLines)+2)[1:-1]

    for line, x, label in zip(labLines, xvals, labels):
        labelLine(line, x, label, align, **kwargs)


# ####################################
# Legend utilities
# ####################################
def legend_yerr(axes, loc='best'):
    """Makes a legend for the ax object 'axes' at the location 'loc',
    but only includes the y error errorbars in the legend icons.

    Parameters
    ----------
    axes : axes.Axes
        The axes on which we want to create a legend with only
        y error bars
    loc : str
        The location of the legend.
    """
    handles, labels = axes.get_legend_handles_labels()

    new_handles = []

    for _, h in enumerate(handles):
        #only need to edit the errorbar legend entries
        if isinstance(h, container.ErrorbarContainer):
            new_handles.append(
                container.ErrorbarContainer(
                    h.lines, has_xerr=False, has_yerr=True
                    )
                )
        else:
            new_handles.append(h)

    axes.legend(new_handles, labels, loc=loc)


def legend_darklight(axes, darklabel='Monte Carlo',
                     lightlabel='Analytic',
                     errtype=None, twosigma=False,
                     extralabel=None):
    """Makes a legend for the ax object 'axes' at the location
    'loc', which indicates dark solid lines as 'MC' (or darklabel)
    and light dotted lines as 'Analytic' (or lightlabel).

    Parameters
    ----------
    axes : axes.Axes
        Description of parameter `axes`.
    darklabel : str
        Label of the dark objects in the legend.
    lightlabel : str
        Label of the light objects in the legend.
    errtype : str
        Specifies the plot type of the `dark` data,
        and thus the corresponding legend icon.
        None (default): legend for line plot.
        yerr: legend for a plot with ecolor
        modstyle: legend for a modstyle errorbar plot
    twosigma : bool
        Determines whether to include lighter, two sigma
        error bars in the legend.

    Returns
    -------
    type
        Description of returned object.

    """
    if errtype is None:
        custom_lines = [Line2D([0], [0], **style_solid,
                               color=compcolors[(-1, 'dark')]),
                        Line2D([0], [0], **style_dashed,
                               color=compcolors[(-1, 'light')])]

        axes.legend(custom_lines, [darklabel, lightlabel])

    elif errtype == 'yerr':
        _, xmax = axes.get_xlim()
        axes.errorbar(xmax*50., 0, yerr=1., **style_yerr,
                      color='black', ecolor='gray', label=darklabel)

        if twosigma:
            axes.errorbar(xmax*50., 0, yerr=2., **style_yerr,
                          color='black', ecolor='lightgray',
                          label=darklabel)
        if extralabel is not None:
            axes.errorbar(xmax*50., 0, yerr=1., **style_yerr_ps,
                          color=compcolors[(-1, 'medium')],
                          ecolor=compcolors[(-1, 'light')],
                          label=extralabel)
            if twosigma:
                axes.errorbar(xmax*50., 0, yerr=2., **style_yerr_ps,
                              color=compcolors[(-1, 'medium')],
                              ecolor=compcolors[(-1, 'light')],
                              label=extralabel)

        handles, _ = axes.get_legend_handles_labels()

        if twosigma:
            l = 0
            # Setting up containers for both errorbars
            if extralabel is not None:
                twosig_extra = container.ErrorbarContainer(
                    handles[-1].lines, has_xerr=False,
                    has_yerr=True)
                onesig_extra = container.ErrorbarContainer(
                    handles[-2].lines, has_xerr=False,
                    has_yerr=True)
                l = -2
            twosig = container.ErrorbarContainer(
                handles[l-1].lines, has_xerr=False,
                has_yerr=True)
            # Setting up containers for both errorbars
            onesig = container.ErrorbarContainer(
                handles[l-2].lines, has_xerr=False,
                has_yerr=True)
            if extralabel is not None:
                custom_handles = [(twosig, onesig),
                                  (twosig_extra, onesig_extra)]
            else:
                custom_handles = [(twosig, onesig)]
        else:
            l = 0
            # Setting up containers for both errorbars
            if extralabel is not None:
                onesig_extra = container.ErrorbarContainer(
                    handles[-1].lines, has_xerr=False,
                    has_yerr=True)
                l=-1
            onesig = container.ErrorbarContainer(
                handles[l-1].lines, has_xerr=False,
                has_yerr=True)
            if extralabel is not None:
                custom_handles = [onesig, onesig_extra]
            else:
                custom_handles = [onesig]

        custom_handles.append(
            Line2D([0], [0], **style_dashed, color=compcolors[(-1, 'light')])
            )

        if twosigma:
            if extralabel is not None:
                axes.legend(custom_handles, [darklabel, extralabel, lightlabel],
                            handler_map={
                                onesig: HandlerErrorbar(xerr_size=0.37),
                                twosig: HandlerErrorbar(xerr_size=0.65),
                                onesig_extra: HandlerErrorbar(xerr_size=0.37),
                                twosig_extra: HandlerErrorbar(xerr_size=0.65)}
                            )
            else:
                axes.legend(custom_handles, [darklabel, lightlabel],
                            handler_map={
                                onesig: HandlerErrorbar(xerr_size=0.37),
                                twosig: HandlerErrorbar(xerr_size=0.65)}
                            )
        else:
            if extralabel is not None:
                axes.legend(custom_handles, [darklabel, extralabel, lightlabel])
            else:
                axes.legend(custom_handles, [darklabel, lightlabel])

    elif errtype == 'modstyle':
        _, xmax = axes.get_xlim()
        axes.errorbar(xmax*50., 0, yerr=1., **modstyle,
                       color=compcolors[(-1, 'dark')],
                       label=darklabel)
        if extralabel is not None:
            axes.errorbar(xmax*50., 0, yerr=1., **modstyle_ps,
                          color=compcolors[(-1, 'medium')],
                          label=extralabel)

        handles, _ = axes.get_legend_handles_labels()

        # Setting up containers for errorbars
        l = 0
        if extralabel is not None:
            onesig_extra = container.ErrorbarContainer(
                handles[-1].lines, has_xerr=True,
                has_yerr=True)
            l=-1
        onesig = container.ErrorbarContainer(handles[l-1].lines,
                                             has_xerr=True,
                                             has_yerr=True)
        if extralabel is not None:
            custom_handles = [onesig, onesig_extra]
        else:
            custom_handles = [onesig]

        custom_handles.append(
            Line2D([0], [0], **style_dashed, color=compcolors[(-1, 'light')])
            )

        if extralabel is not None:
            axes.legend(custom_handles, [darklabel, extralabel, lightlabel],
                        prop={'size': 15}, loc='upper left')
        else:
            axes.legend(custom_handles, [darklabel, lightlabel])
