#!/usr/bin/python3

"""ewoc_plot_utils.py: A set of utilites for reading and plotting
EWOC information from text files."""
__author__ = "Samuel Alipour-fard, Ian Moult, Wouter Waalewijn"
__credits__ = ["Samuel Alipour-fard",
               "Ian Moult",
               "Wouter Waalewijn"]
__license__ = "MIT"
__version__ = "0.3"
__maintainer__ = "Samuel Alipour-fard"
__email__ = "samuelaf@mit.edu"
__status__ = "Production"

# ---------------------------------
# Imports
# ---------------------------------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from cycler import cycler

import inspect
from functools import reduce

from plot_utils import aestheticfig, modstyle

# =====================================
# Utilities for reading EWOC data  
# =====================================
def ewoc_folder(argv):
    """
    Returns a folder used to store subjet pair information given
    a set of command line inputs

    Parameters
    ----------
        argv :  A list containing the command line inputs.

    Returns
    -------
        string : The relevant folder.

    """
    assert len(argv) > 6, "Invalid number of command line parameters.\n"\
        + "    (found " + str(len(argv)) + ", expected 6-8)\n"
    n_events    = argv[1]
    qcd_level   = argv[2]
    process_str = argv[3]
    # no ptmin/max
    jet_alg_int = argv[4]
    jet_rad     = argv[5]

    return "../output/" + process_str + "_" + qcd_level\
        + f"/jetR{float(jet_rad):.1f}/".replace(".", "-")\
        + alg_to_string(jet_alg_int, latex=False) + "jet/"


def ewoc_file_label(argv):
    """
    Returns a file used to store subjet pair information given
    a set of command line inputs

    Parameters
    ----------
        argv :  A list containing the command line inputs.

    Returns
    -------
        string : The relevant file.
    """
    assert len(argv) == 8, "Invalid number of command line parameters.\n"\
        + "    (found " + str(len(argv)) + ", expected 9)\n"
    n_events    = argv[1]
    qcd_level   = argv[2]
    process_str = argv[3]
    # no ptmin/max
    jet_alg_int = argv[4]
    jet_rad     = argv[5]
    sub_alg_int = argv[6]
    sub_rad     = argv[7]

    filename = f"subR{float(sub_rad):.2f}".replace(".", "-") + "_"\
        + alg_to_string(sub_alg_int, latex=False) + "sub"\
        + "_" + str(n_events) + "evts"

    return ewoc_folder(argv) + filename


def get_hist_dict(argv):
    file_prefix = ewoc_file_label(argv)
    data_dir = file_prefix + '.txt'
    hist_dir = file_prefix + '.npz'

    try:
        # Loading data if it exists
        loaded_dict = np.load(hist_dir, allow_pickle=True)
        hist_dict = {}
        for key in loaded_dict.keys():
            hist_dict[key] = loaded_dict[key][()]
        del loaded_dict

    except FileNotFoundError as e:
        # Otherwise, creating the data
        print(e, "Processing data now.")
        data_dict = ewoc_text_to_dict(data_dir)
        hist_dict = ewoc_dict_to_hists(data_dict) 
        np.savez(hist_dir, **hist_dict)

    return hist_dict


def alg_to_string(jet_alg_int, latex=True):
    """
    Returns a string naming the jet algorithm with the
    given index.

    Parameters
    ----------
        jet_alg_int : int
            Integer indicating the FastJet index associated with the jet
            algorithm.

    Returns
    -------
        str
            A LaTeX-compatible string naming the algorithm.
    """
    if jet_alg_int in [0, "0"]:
        if latex:
            return r'$k_T$'
        else:
            return 'kt'
    elif jet_alg_int in [1, "1"]:
        if latex:
            return r'C/A'
        else:
            return 'ca'
    elif jet_alg_int == [2, "2"]:
        if latex:
            return r'anti-$k_T$'
        else:
            return 'akt'
    else:
        raise AssertionError("Invalid jet algorithm index " +
                         str(jet_alg_int))


# =====================================
# Pairwise subjet observables (collinear limit) 
# =====================================

def pair_costheta(Etot, z1, z2, costheta):
    return costheta


def pair_z(Etot, z1, z2, costheta):
    return (1.-costheta)/2.


def pair_formtime(Etot, z1, z2, costheta):
    try:
        return Etot * max(z1, z2) / pair_mass(Etot, z1, z2, costheta)
    except ZeroDivisionError:
        return np.nan



def pair_mass(Etot, z1, z2, costheta):
    return 2. * Etot**2. * z1 * z2 * (1.-costheta)


# =====================================
# Formatting Data 
# =====================================

def ewoc_text_to_dict(filename, pair_obs=None):
    """
    Reads in a text file in the form provided by pythia_ewocs(.cc).
    Returns a dict whose highest keys are 'info' and a list of subjet radii,
    and whose 2nd level keys are properties -- 'weights' (z1*z2),
    'costheta', 'mass', and 'formtime'. Roughly,
        ```data_dict[property] = list```

    Parameters
    ----------
        filename : str
            Name of an ewoc formatted file.

    Returns
    -------
        data_dict : dict
            A dict containing organized information from the ewoc file. 
    """
    # ---------------------------------
    # Setup 
    # ---------------------------------
    # Observables
    if pair_obs is None:
        pair_obs = [pair_costheta, pair_formtime, pair_mass]
    if not isinstance(pair_obs, list):
        # Making sure observables are stored in a list
        pair_obs = [pair_obs]
    # storing observable names
    obsnames = [f'{obsname=}'.split(' ')[1][len("pair_"):]
                for obsname in pair_obs]

    # Event headers/information: Event, Jet, or Subjet Pair
    event_headers = ['E', 'J', 'SP']
    i_event = 0
    E_jet = -1

    # Preparing to store data from text file
    data_dict = {obsname: [] for obsname in obsnames}
    data_dict['weights'] = []
    data_dict['observables'] = obsnames

    # ---------------------------------
    # Reading file
    # ---------------------------------
    with open(filename, "r") as file:
        for line in file:
            # Read and format the line
            line = line.rstrip("\n")
            info = line.split(" ")

            # Additional information:
            if info[0] == "#" or '' in info:
                # Allowing commented or empty lines
                continue
            elif info[0] not in event_headers:
                # Storing additional information from file header
                data_dict[info[0]] = info[-1]
                continue

            # Event information
            if info[0] == "E":
                if i_event%1000 == 0 and i_event > 0:
                    # Better logging
                    print(f"Considered {i_event} events")
                i_event+=1

            if info[0] == "J":
                E_jet = float(info[1])

            if info[0] == "SP":
                # Convert strings to floats
                z1, z2, costheta = (float(x) for x in info[1:])

                # Store relevant EWOC information
                data_dict['weights'].append(z1 * z2)
                for obsname, obs in zip(obsnames, pair_obs):
                    data_dict[obsname].append(obs(E_jet, z1, z2, costheta))
                
    # After the loop, put weights and observables into numpy arrays
    for obsname in ['weights', *obsnames]:
        data_dict[obsname] = np.array(data_dict[obsname])

    return data_dict


def get_hist_edges_centers(obs_vals, weights, nbins, binspace):
    """
    Takes in a set of values and an associated binning space.
    Returns (histogram heights, bin edges, bin centers)
    associated with those values.

    Parameters
    ----------
        obs_vals :  The set of values to histogram.
        nbins :     The number of desired bins.
        binspace :  The desired bin space ('lin'/'linear' or 'log').

    Returns
    -------
        (hist, bin_edges, bin_centers) :
            Three lists describing the properties of the histogram.

    """
    lin_labels, log_labels = ['linear', 'lin'], ['log']

    o_vals = obs_vals if binspace in lin_labels\
             else np.log(obs_vals) if binspace in log_labels\
             else None

    # Making and storing weighted and normalized histogram
    hist, bin_edges = np.histogram(o_vals,
            np.linspace(np.min(o_vals), np.max(o_vals), nbins+1),
            weights=weights, density=True)

    if binspace in ['linear', 'lin']:
        bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.

    elif binspace in ['log']:
        # Getting logarithmic hists for appropriate observables
        bin_centers = np.exp((bin_edges[:-1] +
                              bin_edges[1:]) / 2.)
        bin_edges = np.exp(bin_edges)

    return (hist, bin_edges, bin_centers)


def ewoc_dict_to_hists(data_dict, hist_dict=None,
                       nbins=500, binspaces=['linear', 'log']):
    """
    Takes in the location of an ewoc formatted file.
    Returns a dict with (roughly)
        ```dict[subjet_radius][observable][lin or log] =
            (histogram values, hist bin edges)
        ```
    for costheta, z = (1-costheta)/2, and mass ewocs.

    Parameters
    ----------
        filename : str
            Name of the ewoc formatted file.
        nbins : int
            Number of bins (not bin edges!) for the hists.

    Returns
    -------
    Dict of histograms containing easy-to-plot ewoc information.
    """
    # Preparing to store histograms for this info
    if hist_dict is None:
        hist_dict = {'observables': data_dict['observables']}
    
    # Adding hists/bin edges to a dict for each subjet radius/observable
    for key, data in data_dict.items():
        if not isinstance(data, np.ndarray):
            hist_dict[key] = data
        elif key in data_dict['observables']:
            if key == 'costheta':
                # cos(theta) comes along with z=(1-cos)/2:
                hist_dict[key] = {space:
                                  get_hist_edges_centers(data_dict[key],
                                                         data_dict['weights'],
                                                         nbins, space)
                                  for space in binspaces
                                  if space in ['linear', 'lin']}

                # z = (1-cos(theta))/2
                hist_dict['z'] = {}
                hist_dict['observables'].append('z')
                hist_dict['z'] = {space:
                                  get_hist_edges_centers((1.-data_dict[key])/2.,
                                                         data_dict['weights'],
                                                         nbins, space)
                                  for space in binspaces}
                continue
            hist_dict[key] = {space:
                              get_hist_edges_centers(data_dict[key],
                                                     data_dict['weights'],
                                                     nbins, space)
                              for space in binspaces}

    return hist_dict


# =====================================
# Plotting Utilites 
# =====================================

# ---------------------------------
# Plotting utils
# ---------------------------------

lims = {'costheta': {('linear', 'linear'): 
                        #((-1.0, 1.0), (0, 30)),
                        ((0.6, 1.0), (0, 10)),
                      ('linear', 'log'): 
                        ((.65, 1.0), (.1, 50)),
                      ('log', 'linear'):
                        ((0.0, 0.5), (0, 30)),
                      ('log', 'log'):
                        ((.65, 1.0), (.1, 50))
                      },
         'z'       : {('linear', 'linear'): 
                        ((0, .25), (0, 20)),
                      ('linear', 'log'): 
                        ((0, .25), (.1, 50)),
                      ('log', 'linear'):
                        ((1e-6, 1e0), (0, .4)),
                      ('log', 'log'):
                        ((1e-6, .2), (.1, 50))
                      },
         'mass'   : {('linear', 'linear'): 
                        ((0, 5e5), (0, 1e-4)),
                      ('linear', 'log'): 
                        ((0, 150), (1e-4, .12)),
                      ('log', 'linear'):
                        ((1e-2, 5e5), (0, .4)),
                      ('log', 'log'):
                        ((1.0, 1e3), (1e-4, .7))
                      },
        'formtime': {('linear', 'linear'): 
                        ((0, 1e7), (0, 2e0)),
                      ('linear', 'log'): 
                        ((0, 150), (1e-4, .12)),
                      ('log', 'linear'):
                        ((1e-4, 1e5), (0, 1e0)),
                      ('log', 'log'):
                        ((1.0, 1e3), (1e-4, .7))
                      }
        }


# Colormap for smooth
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


# =====================================
# Old 
# =====================================

def obs_label(observable):
    """Returns a label for plotting the given observable
    in LaTeX form.
    """
    if observable == 'costheta':
        return r'$\cos(\theta)$ '
    elif observable == 'z':
        return r'$z = (1 - \cos\theta)/2$ '
    elif observable == 'mass':
        return r'Mass (GeV)'
    elif observable == 'formtime':
        return r'Formation Time (GeV$^{-1}$)'


def hist_title(observable, jet_radius, subjet_alg_int):
    """Title for a histogram plotting the EWOC for the given
    subjet pair property, jet radius, and subjet algorithm.
    """
    title = ''
    """
    if observable == 'costheta':
        title += r'$\cos(\theta)$ '
    elif observable == 'z':
        title += r'$z = (1 - \cos\theta)/2$ '
    elif observable == 'mass':
        title += r'Mass '
    elif observable == 'formtime':
        title += r'Formation Time '
    """

    # radius info
    title += r'EWOC, $R_{\rm jet}$=' + f'{float(jet_radius):.1f}'

    # Add subjet alg info
    title += f', {alg_to_string(subjet_alg_int)} Subjets'
            

    return title


# Utils for QCD predictions
M_Z = 91.19 # GeV

# Group theory factors
CF = 8./6.
CA = 3.
TF = 1./2.
N_F = 5.

# QCD Beta function 
beta_0 = (33 - 2*N_F)/(12 * np.pi)


def alpha_s(mu):
    """1-loop strong force coupling. Argument in GeV."""
    return .12/(1 + 2*.12*beta_0*np.log(mu/M_Z))


def eec_nnlo(z, mu):
    """NNLO EEC."""
    a_s = alpha_s(mu) / (4 * np.pi)
    # Not quite right yet, need to run at higher acccuracy
    
    z = z.astype(float)

    nlo_piece  = -11.5333*np.log(z) + 81.4809
    nnlo_piece = 45.1489*np.log(z)**2. - 1037.73*np.log(z) + 2871.36

    return (2.*a_s + a_s**2.*nlo_piece + a_s**3.*nnlo_piece)/(z*(1-z))

# ---------------------------------
# Plotting functions
# ---------------------------------

def plot_eec_analytic(ax, observable, binspace, acc='nnlo'):
    """Plots an analytic EEC to the given accuracy on the given axes."""
    assert acc == 'nnlo', "Invalid accuracy " + str(acc)

    # rough, for now
    Q = 2000

    # Plot the full analytic EEC
    if observable == 'zs':
        xs = np.concatenate((np.linspace(1e-8, 1, 250),
                             np.logspace(-8, 0, 250)))
        xs = np.sort(xs)
        ax.plot(xs, eec_nnlo(xs, Q), 'k--',
                linewidth=2, zorder=4)
    if observable == 'costhetas':
        xs = np.linspace(-1, 1, 500)
        ax.plot(xs, eec_nnlo((1-xs)/2, Q), 'k--',
                linewidth=2, zorder=4)
    return


def plot_hist_dict(hist_dicts, obsname,
                   binspace=None, ax=None,
                   xlim='default',
                   ylim='default',
                   show_analytic=True):
    # Showing analytic EWOC for only specific observables
    show_analytic = obsname in ['costheta', 'z'] and show_analytic

    if ax is None:
        # Set up figure, axis limits
        _, ax = aestheticfig(xlabel=obs_label(obsname),
                             ylabel='EWOC',
                             title=None,
                             showdate=False,
                             xlim=xlim if xlim != 'default' else
                                  lims[obsname][(binspace, 'linear')][0],
                             ylim=ylim if ylim != 'default' else
                                  lims[obsname][(binspace, 'linear')][1],
                             ratio_plot=(obsname in ['costheta',
                                                     'z']),
                             labeltext='EWOC')

    if not isinstance(hist_dicts, list):
        hist_dicts = [hist_dicts]

    # Looping over subjet radii
    for hist_dict in hist_dicts:
        # Set up x errorbar
        ys, x_edges, xs = hist_dict[obsname][binspace]
        ax[0].errorbar(xs, ys, yerr=np.sqrt(ys/len(ys)),
                       xerr=(xs - x_edges[:-1], x_edges[1:] - xs),
                       **modstyle)
        ax[0].set_title(hist_title(obsname, hist_dict['jet_rad'],
                                   hist_dict['sub_alg']))

        if show_analytic:
            analytic = eec_nnlo((1.-xs)/2., 2000)\
                if obsname == 'costheta'\
                else eec_nnlo(xs, 2000)

            ax[1].errorbar(xs, ys/analytic,
                           yerr=np.sqrt(ys/len(ys))/analytic,
                           xerr=(xs - x_edges[:-1], x_edges[1:] - xs),
                           **modstyle)

    if show_analytic:
        ax[0].plot(xs, analytic, 'k')
        ax[1].plot(xs, np.ones(len(xs)), 'k')
        ax[1].set_ylim((1e-1, 1e1))

    #plot_utils.stamp(text)
    if binspace == 'log':
        [a.set_xscale('log') for a in ax]

    return


# def plot_single_subR(argv):
#     file_prefix = ewoc_file_label(argv)

#     data_dir = file_prefix + '.txt'
#     hist_dir = file_prefix + '.npz'

#     try:
#         # Loading data if it exists
#         loaded_dict = np.load(hist_dir, allow_pickle=True)
#         hist_dict = {}
#         for key in loaded_dict.keys():
#             hist_dict[key] = loaded_dict[key][()]
#         del loaded_dict

#     except FileNotFoundError as e:
#         # Otherwise, creating the data
#         print(e.message())
#         data_dict = ewoc_text_to_dict(data_dir)
#         hist_dict = ewoc_dict_to_hists(data_dict) 
#         np.savez(hist_dir, **hist_dict)

#     # Plotting
#     for obsname in obsnames:
#         # yspace -- use itertools
#         for binspace in ['linear', 'log']:
#             # No logarithmic bins for costheta
#             if binspace == 'log' and obsname == 'costheta':
#                 continue
#             plot_hist_dict(hist_dict)

#     return


def plot_EWOC_by_sub_rads(argv):
    # Getting list of radii from command line argument
    subjet_rads = argv[7].lstrip('[').rstrip(']').split(',')
    subjet_rads = [float(R) for R in subjet_rads]
    colors, cbar_data = get_colors_colorbar(subjet_rads)

    # Getting observable dictionaries
    hist_dicts = []
    for sub_rad in subjet_rads:
        hist_dicts.append(get_hist_dict([*argv[:7], sub_rad]))

    # Finding the shared observables in each dict
    obsnames = list(reduce(set.intersection,
                           [set(hist_dict['observables'])
                            for hist_dict in hist_dicts])
                   )

    # Plotting
    for obsname in obsnames:
        # yspace -- use itertools
        for binspace in ['linear', 'log']:
            # No logarithmic bins for costheta
            if binspace == 'log' and obsname == 'costheta':
                continue

            # Set up figure, axis limits
            fig, ax = aestheticfig(xlabel=obs_label(obsname),
                                   ylabel='EWOC',
                                   title=None,
                                   showdate=False,
                                   xlim=lims[obsname][(binspace, 'linear')][0],
                                   ylim=lims[obsname][(binspace, 'linear')][1],
                                   ratio_plot=(obsname in ['costheta',
                                                           'z']),
                                   labeltext='EWOC')
            # Setting up color cycle
            [a.set_prop_cycle((cycler(color=colors))) for a in ax]

            # Plotting observable EWOC for each subjet radius
            plot_hist_dict(hist_dicts, obsname,
                           binspace=binspace,
                           ax=ax)

            # Adding subjet radius color bar
            cbar = fig.colorbar(cbar_data, ax=np.array(ax).ravel().tolist())
            cbar.ax.set_ylabel(r'$R_{\rm sub}$')

    plt.show()

    return
    """Plots a set of EWOCs associated with a given
    dictionary.
    """
    """
    # Older example
    if obsnames is None:
        obsnames = hist_dict['observables']

    # Grabbing subjet radii (numbers)
    #subjet_rads = [rad for rad in hist_dict.keys()
    #               if not isinstance(rad, str)]
    # Setting up colors for plotting
    # colors, cbar_data = get_colors_colorbar(subjet_rads)

    # Plotting ewocs for each property in log and lin space
    for propname in props:
        for binspace in ['lin', 'log']:
            xscale = 'linear' if binspace=='lin' else binspace
            # Plotting linear and logarithmic observables
            if propname == 'costhetas' and binspace == 'log':
                continue
            for yscale in ['linear', 'log']:
                # Using both linear and logarithmic y axes
                fig, ax = plt.subplots()
                for ic, sub_rad in enumerate(subjet_rads):
                    # Looping over subjet radii
                    ys, edges, xs = hist_dict[sub_rad][propname][binspace] 
                    ax.plot(xs, ys, label=r'$R_{\rm sub}$='+str(sub_rad),
                            color=colors[ic])

                    # Setting axis limits
                    if xlim == 'default':
                        ax.set_xlim(lims[propname][(xscale, yscale)][0])
                    elif xlim is not None:
                        ax.set_xlim(xlim)
                    elif xscale == 'log':
                        ax.set_xlim(np.nanmin([e for e in edges if e>0]),
                                    min(np.max(edges), 1e4))
                    if ylim == 'default':
                        ax.set_ylim(lims[propname][(xscale, yscale)][1])
                    elif ylim is not None:
                        ax.set_ylim(ylim)
                    elif yscale == 'log':
                        ax.set_ylim(np.nanmin([y for y in ys if y>0]),
                                    np.max(ys))
                if show_analytic and (propname == 'costhetas' or
                                      propname == 'zs'):
                    plot_eec_analytic(fig, ax, propname, binspace)
                # Other plot setup
                try:
                    ax.set_xscale(xscale)
                    ax.set_yscale(yscale)
                except ValueError as e:
                    # Data with no positive values cannot be log scaled
                    print(f"Error message for {propname},"
                          + f"xscale={xscale}, yscale={yscale}:\n\n")
                    print(e)
                    continue
                #ax.legend()
                cbar = fig.colorbar(cbar_data)
                cbar.ax.set_ylabel(r'$R_{\rm sub}$')

                ax.set_title(hist_title(propname, hist_dict['radius'],
                                          hist_dict['subjet algorithm']))

                fig.savefig('figs/' + hist_dict['process'] + '_'
                            + hist_dict['level'] + '/'
                            # e.g. figs/qcd_parton
                            + propname+"_ewoc_sub"
                            + hist_dict['subjet algorithm']
                            + '_' + xscale + 'x'
                            + '_' + yscale + 'y'
                            + "_R" +
                            str("{:.1f}".format(float(hist_dict['radius'])))
                            + "_" + str(hist_dict['nEvents'])
                            + "_" + hist_dict['level'] + "_"
                            + hist_dict['process'] + 'evts.pdf',
                            bbox_inches='tight')

                # plt.show()
                
    return
    """
