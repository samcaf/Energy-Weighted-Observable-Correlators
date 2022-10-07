#!/usr/bin/python3

"""ewoc_utils.py:
A set of utilites for reading and plotting
EWOC information from text files.
"""

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
import itertools
from functools import reduce

import os.path

from plot_utils import aestheticfig, aesthetic_N_by_M, stamp
from plot_utils import style_solid, style_dashed, style_scatter, modstyle

from ewoc_cmdln import default_args, arg_to_str_dict, arg_to_list

from matplotlib import rc
rc('text', usetex=True)


# ---------------------------------
# Global Flags
# ---------------------------------
plot_errorbars = False
plot_scatter = True  # Only applied if plot_errorbars = False


# =====================================
# Utilities for reading EWOC data  
# =====================================
def ewoc_folder(n_events, qcd_level, process_str,
        jet_alg_int, jet_rad, *args, **kwargs):
    """
    Returns a folder used to store subjet pair information given
    a set of command line inputs

    Parameters
    ----------
        Command line arguments; see ewoc_cmdln.py

    Returns
    -------
        string : The relevant folder.

    """
    return "../output/" + process_str + "_" + qcd_level\
        + f"/jetR{float(jet_rad):.1f}/".replace(".", "-")\
        + alg_to_string(jet_alg_int, latex=False) + "jet/"


def ewoc_file_label(n_events, qcd_level, process_str,
        jet_alg_int, jet_rad, sub_alg_int, sub_rad,
        **kwargs):
    """
    Returns a file used to store subjet pair information given
    a set of command line inputs

    Parameters
    ----------
        Command line arguments; see ewoc_cmdln.py

    Returns
    -------
        string : The relevant file.
    """
    folder = ewoc_folder(n_events, qcd_level, process_str,
                jet_alg_int, jet_rad, sub_alg_int, sub_rad,
                **kwargs)
    filename = f"subR{float(sub_rad):.2f}".replace(".", "-") + "_"\
        + alg_to_string(sub_alg_int, latex=False) + "sub"\
        + "_" + str(n_events) + "evts"

    for key in arg_to_str_dict.keys():
        if key in kwargs.keys():
            if kwargs[key] != default_args[key] and kwargs[key] is not None: 
                filename += '_' + arg_to_str_dict[key]\
                         + str(kwargs[key]).replace(".", "-")

    return folder + filename


def save_hist_dict(overwrite=False, print_every_n=1000,
                  **kwargs):
    file_prefix = ewoc_file_label(**kwargs)
    data_dir = file_prefix + '.txt'
    # DEBUG: Using to generate new files when attempting
    # to normalize ewocs by hand
    hist_dir = file_prefix + '_normed.npz'
    #hist_dir = file_prefix + '.npz'

    if os.path.isfile(hist_dir):
        print(f"File found at {hist_dir}.")
        if not overwrite:
            print("Not overwriting existing file.\n")
            return
        else:
            print("Overwriting existing file.\n")
    else:
        print(f"No file found at {hist_dir}.")

    # Creating the histogram data
    print(f"Processing data from {data_dir}.\n")
    data_dict = ewoc_text_to_dict(data_dir, print_every_n=print_every_n)
    hist_dict = ewoc_dict_to_hists(data_dict) 
    del data_dict

    print(f"Saving histogram data to {hist_dir}"\
          + " and deleting stored dictionary.")
    np.savez(hist_dir, **hist_dict)
    del hist_dict

    return


def get_hist_dict(load=True, save=True, print_every_n=1000,
                  **kwargs):
    file_prefix = ewoc_file_label(**kwargs)
    data_dir = file_prefix + '.txt'
    # DEBUG: Using to generate new files when attempting
    # to normalize ewocs by hand
    hist_dir = file_prefix + '_normed.npz'
    #hist_dir = file_prefix + '.npz'

    if load: 
        try:
            # Loading data if it exists
            print(f"Loading histogram data from {hist_dir}.")
            loaded_dict = np.load(hist_dir, allow_pickle=True)
            hist_dict = {}

            for key in loaded_dict.keys():
                # Respect weird syntax from dicts saved w/```np.savez```
                hist_dict[key] = loaded_dict[key][()]
            del loaded_dict

            return hist_dict

        except FileNotFoundError as e:
            print(f"Unable to find {hist_dir}:\n   ", e)

    # Otherwise, creating the histogram data
    print(f"Processing data from {data_dir}.")
    data_dict = ewoc_text_to_dict(data_dir, print_every_n=print_every_n)
    hist_dict = ewoc_dict_to_hists(data_dict) 
    del data_dict

    if save:
        print(f"Saving histogram data to {hist_dir}.")
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
    elif jet_alg_int in [2, "2"]:
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

def obs_label(observable):
    """Returns a label for plotting the given observable
    in LaTeX form.
    """
    if observable == 'costheta':
        return r'$\cos(\theta)$'
    elif observable == 'z':
        return r'$z = (1 - \cos\theta)/2$'
    elif observable == 'mass':
        return r'Mass (GeV)'
    elif observable == 'formtime':
        return r'Formation Time (GeV$^{-1}$)'

def obs_title(observable):
    """Returns a title for plotting the given observable
    in LaTeX form.
    """
    if observable == 'costheta':
        return r'$\cos(\theta)$ Subjet EWOC'
    elif observable == 'z':
        return r'$z$ Subjet EWOC'
    elif observable == 'mass':
        return r'Mass Subjet EWOC'
    elif observable == 'formtime':
        return r'Formation Time Subjet EWOC'


def pair_costheta(Etot, z1, z2, costheta):
    return costheta

def pair_z(Etot, z1, z2, costheta):
    return (1.-costheta)/2.

def pair_m2(Etot, z1, z2, costheta):
    return 2. * Etot**2. * z1 * z2 * (1.-costheta)

def pair_mass(Etot, z1, z2, costheta):
    return np.sqrt(pair_m2(Etot, z1, z2, costheta))

def pair_formtime(Etot, z1, z2, costheta):
    return np.float64(Etot * max(z1, z2)) / pair_m2(Etot, z1, z2, costheta)
    # try:
    # except ZeroDivisionError:
    #     return np.nan


# =====================================
# Formatting Data 
# =====================================

def ewoc_text_to_dict(filename, pair_obs=None,
                      print_every_n=1000):
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
        pair_obs : list (of functions)
            Functions of E_jet, z1, z2, and cos(theta) corresponding to
            observables to include in the dict.
        print_every_n : int
            Print a message every time ```print_every_n``` events have
            been considered.

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
    # Storing observable names, assuming they start with "pair_"
    obsnames = [f'{obsname=}'.split(' ')[1][len("pair_"):]
                for obsname in pair_obs]

    # Event headers/information: Event, Jet, or Subjet Pair
    event_headers = ['E', 'J', 'SP']
    i_event = 1
    E_jet = -1

    # Preparing to store data from text file
    data_dict = {obsname: [] for obsname in obsnames}
    data_dict['weights'] = []
    data_dict['observables'] = obsnames

    # ---------------------------------
    # Reading file
    # ---------------------------------
    with open(filename, "r") as file:
        for i, line in enumerate(file):
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

            if info[0] == "J":
                E_jet = float(info[1])

            if info[0] == "SP":
                # Convert strings to floats
                try:
                    z1, z2, costheta = (float(x) for x in info[1:])
                except ValueError as e:
                    print(e)
                    print(f"    line in file : {info}")
                    print(f"    Next line in file: {file[i+1]}")
                    raise ValueError

                # Store relevant EWOC information
                data_dict['weights'].append(z1 * z2)
                for obsname, obs in zip(obsnames, pair_obs):
                    data_dict[obsname].append(obs(E_jet, z1, z2, costheta))

            # Event information
            if info[0] == "E":
                if i_event%print_every_n == 0:
                    # Better logging
                    print(f"Considered {i_event} events")
                i_event+=1
               
    # After the loop, put weights and observables into numpy arrays
    for obsname in ['weights', *obsnames]:
        data_dict[obsname] = np.array(data_dict[obsname])

    return data_dict


def normalized_hist(obs_vals, bin_edges, weights, binspace):
    """
    Returns a histogram for the given which is normalized over
    the full range of those values, but only within the bin_edges
    provided. Allows for unconvenional "logarithmic" normalization.

    Parameters
    ----------
        obs_vals :      The observables to histogram;
        bin_edges :     The edges of the histogram bins;
        weights :       Weights of the histogram;
        binspace :      Parameter to determine normalization:
                        if 'lin' or 'linear', normalizes in linear
                        space;
                        if 'log', assumes the observable and bins are
                        logarithms of the "true" values, and normalizes
                        in the "true" space.

    Returns
    -------
       A histogram, normalized over the range of obs_vals, for obs_vals
       within bin_edges.
    """
    assert binspace in ['lin', 'linear', 'log'], f"Invalid {binspace=}"
    if binspace in ['lin', 'linear']:
        full_hist, full_bins = np.histogram(obs_vals,
                            np.linspace(np.nanmin(obs_vals), np.nanmax(obs_vals), 2),
                            weights=weights)
        integral = np.sum(hist * (full_bins[1:] - full_bins[:-1]))
        hist = np.histogram(obs_vals, bin_edges, weights=weights/integral)

    elif binspace == 'log':
        full_hist, full_bins = np.histogram(obs_vals,
                            np.linspace(np.nanmin(obs_vals), np.nanmax(obs_vals), 2),
                            weights=weights)

        hist = None
        # DEBUG
        

    return hist
    

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

    # DEBUG
    #hist, _ = normalized_hist(o_vals, bin_edges, weights, binspace)
    hist, _ = np.histogram(o_vals, bin_edges,
                            weights=weights,
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


binedge_dict = {'mass': {'lbin': 1e-5, 'rbin': 2e3},
                'formtime': {'lbin': 5e-4, 'rbin': 1e4}}


def ewoc_dict_to_hists(data_dict, hist_dict=None,
                       nbins=250, binspaces=['linear', 'log']):
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
                                                     nbins, space)  #, **binedge_dict[key])
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
                        ((3e-5, 1e0), (0, .5)),
                      ('log', 'log'):
                        ((1e-6, .2), (.1, 50))
                      },
         'mass'   : {('linear', 'linear'): 
                        ((0, 1e3), (0, 1.7e-2)),
                      ('linear', 'log'): 
                        ((0, 150), (1e-4, .12)),
                      ('log', 'linear'):
                        ((1e-1, 5e3), (0, .9)),
                      ('log', 'log'):
                        ((1.0, 1e3), (1e-4, .7))
                      },
        'formtime': {('linear', 'linear'): 
                        ((0, 1e4), (0, 1e-4)),
                      ('linear', 'log'): 
                        ((0, 150), (1e-4, .12)),
                      ('log', 'linear'):
                        ((1e-3, 1e3), (0, 4e-1)),
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

def stamp_ax(ax, jet_rad, sub_rad,
             obsname, stamp_loc=None, **kwargs):
    """Add a stamp to the given axes."""
    # Location of the stamp
    if stamp_loc is None:
        stamp_loc = (0.5, 0.95)

    # Text for the stamp
    stamp_text = {'textops_update': {'fontsize': 11,
                      'horizontalalignment': 'center'},
               'line_0':
                f"({alg_to_string(kwargs['jet_alg_int'],latex=False).upper()}"\
                +(f" {jet_rad}" if jet_rad is not None else '') + ", "\
                +f"{alg_to_string(kwargs['sub_alg_int'],latex=False).lower()}"\
                +(f" {sub_rad}" if sub_rad is not None else '')\
                +fr") \textbf{{{obs_title(obsname)}}}",
               'line_1':
                f"{kwargs['n_events']} {kwargs['process_str'].upper()} events, "\
                +r"\texttt{Pythia 8.307}",
               'line_2':
                str(kwargs['pt_min']) + r" $< p_T <$ "\
                + str(kwargs['pt_max'])
            }

    # Make the stamp!
    stamp(*stamp_loc, ax, **stamp_text)

# =====================================
# QCD
# =====================================

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
                   ax=None,
                   binspace=None, y_scale='linear', 
                   xlim='default',
                   ylim='default',
                   title=None,
                   show_analytic=True,
                   hist_weight=None,
                   label=None,
                   **kwargs):
    # Showing analytic EWOC for only specific observables
    show_analytic = obsname in ['costheta', 'z'] and show_analytic

    # Extra weighting factors for plotting the given histogram
    if hist_weight is None:
        # Default: no additional weight
        hist_weight = lambda x: np.ones(len(x))

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
                             ratio_plot=(obsname in ['costheta', 'z']))

    if not isinstance(hist_dicts, list):
        hist_dicts = [hist_dicts]

    # Looping over subjet radii
    for hist_dict in hist_dicts:
        # Set up x errorbar
        ys, x_edges, xs = hist_dict[obsname][binspace]
        ys = ys * hist_weight(xs)

        # Plotting
        if plot_errorbars:
            ax[0].errorbar(xs, ys, yerr=np.sqrt(ys/len(ys)),
                       xerr=(xs - x_edges[:-1], x_edges[1:] - xs),
                       label=label, **modstyle)
        elif plot_scatter:
            ax[0].scatter(xs, ys, label=label, **style_scatter)
        else:
            ax[0].plot(xs, ys, label=label, **style_solid)

        if title is not None:
            ax[0].set_title(title if title != 'default' else\
                hist_title(obsname, hist_dict['jet_rad'],
                           hist_dict['sub_alg']))

        if show_analytic:
            analytic = eec_nnlo((1.-xs)/2., 2000)\
                if obsname == 'costheta'\
                else eec_nnlo(xs, 2000)
            analytic = analytic * hist_weight(xs)

            # log-normalized
            if binspace == 'log':
                analytic = analytic * xs

            if plot_errorbars:
                ax[1].errorbar(xs, ys/analytic,
                           yerr=np.sqrt(ys/len(ys))/analytic,
                           xerr=(xs - x_edges[:-1], x_edges[1:] - xs),
                           **modstyle)
            elif plot_scatter:
                ax[1].scatter(xs, ys/analytic, **style_scatter)
            else:
                ax[1].plot(xs, ys/analytic, **style_solid)

    if show_analytic:
        ax[0].plot(xs, analytic, 'k', **style_dashed)
        ax[1].plot(xs, np.ones(len(xs)), 'k', **style_dashed)
        ax[1].set_ylim((5e-2, 2e1))

    if binspace == 'log':
        [a.set_xscale('log') for a in ax]
    if y_scale == 'log':
        ax[0].set_ylim(bottom=1e-8)
        ax[0].set_yscale('log')

    return


def plot_EWOC_by_sub_rads(load=True, print_every_n=1000,
                          hist_weight=None,
                          **kwargs):
    # Typesetting keyword arguments
    arg_to_list(kwargs, 'sub_rad')
    # input_dict['jet_rad'] = float(input_dict['jet_rad'])

        
    # Saving data
    for rsub in kwargs['sub_rad']:
        save_hist_dict(overwrite=False, print_every_n=print_every_n,
                       **dict(kwargs, sub_rad=rsub))
    
    # Loading and processing data all given subjet radii
    # (accessed via ```kwargs['sub_rad']```)
    hist_dicts = []

    for rsub in kwargs['sub_rad']:
        hist_dicts.append(get_hist_dict(**dict(kwargs, sub_rad=rsub),
                                        load=load,
                                        print_every_n=print_every_n))

    # Finding the shared observables in each dict
    obsnames = list(reduce(set.intersection,
                           [set(hist_dict['observables'])
                            for hist_dict in hist_dicts])
                   )

    # Setting up colorbar
    colors, cbar_data = get_colors_colorbar(kwargs['sub_rad'])

    # Plotting
    for obsname in obsnames:
        for binspace, y_scale in zip(['linear', 'log'],
                                     ['linear', 'log']):
            # No logarithmic bins for costheta
            if binspace == 'log' and obsname == 'costheta':
                continue
            # Don't plot lin-log
            if binspace in ['linear', 'lin'] and y_scale == 'log':
                continue

            # Set up figure, axis limits
            fig, ax = aestheticfig(xlabel=obs_label(obsname),
                                   ylabel='EWOC',
                                   title=None,
                                   showdate=False,
                                   xlim=lims[obsname][(binspace, 'linear')][0],
                                   ylim=lims[obsname][(binspace, 'linear')][1],
                                   ratio_plot=(obsname in ['costheta', 'z']))

            # Setting up color cycle
            [a.set_prop_cycle((cycler(color=colors))) for a in ax]

            # Plotting observable EWOC for each subjet radius
            plot_hist_dict(hist_dicts, obsname, ax=ax,
                           binspace=binspace,
                           y_scale=y_scale,
                           hist_weight=hist_weight,
                           **kwargs)
            if kwargs['process_str'] == 'w':
                if obsname == 'mass':
                    ax[0].plot([80.4, 80.4], [-1, 100], 'r--')
                if obsname == 'formtime':
                    ax[0].plot([1/2.141, 1/2.141], [-1, 100], 'r--')

            # Adding subjet radius color bar
            cbar = fig.colorbar(cbar_data, ax=np.array(ax).ravel().tolist())
            cbar.ax.set_ylabel(r'$R_{\rm sub}$')

            # Stamp text
            stamp_ax(ax[0], jet_rad=kwargs['jet_rad'],
                     sub_rad=None, obsname=obsname,
                     **{key: val for key, val in kwargs.items()
                            if key not in ['jet_rad', 'sub_rad']})

    plt.show()

    return


def plot_EWOC_pvh_by_rads(load=True, print_every_n=1000,
                          hist_weight=None,
                          **kwargs):
    # Turning the keyword argument for (sub)jet radii into a list
    for rad_type in ['sub_rad', 'jet_rad']:
        arg_to_list(kwargs, rad_type)
        assert len(kwargs[rad_type]) == 2, "Expected two radii"\
            + f"of type {rad_type} (received {len(kwargs[rad_type])})."

    # Getting observable dictionaries for all given subjet radii
    # (accessed via ```kwargs['sub_rad']```)
    hist_dicts_parton = []
    hist_dicts_hadron = []
    
    # Saving data
    for (Rjet, rsub, level) in itertools.product(
                                kwargs['jet_rad'],
                                kwargs['sub_rad'],
                                ['parton', 'hadron']):
        save_hist_dict(overwrite=False, print_every_n=print_every_n,
                       **dict(kwargs, qcd_level=level,
                              jet_rad=Rjet, sub_rad=rsub))
    
    # Loading and processing data
    for (Rjet, rsub) in itertools.product(kwargs['jet_rad'],
                                          kwargs['sub_rad']):
        print(f"\n{(Rjet, rsub) = }\n# =====================\n")
        print("Partons:")
        hist_dicts_parton.append(get_hist_dict(
                            **dict(kwargs, qcd_level='parton',
                                   jet_rad=Rjet, sub_rad=rsub),
                            load=load,
                            print_every_n=print_every_n))
        print("Hadrons:")
        hist_dicts_hadron.append(get_hist_dict(
                            **dict(kwargs, qcd_level='hadron',
                                   jet_rad=Rjet, sub_rad=rsub),
                            load=load,
                            print_every_n=print_every_n))
        print()

    # Finding the shared observables in each dict
    obsnames = list(reduce(set.intersection,
                           [set(hist_dict['observables'])
                            for hist_dict in
                            np.concatenate((hist_dicts_parton,
                                            hist_dicts_hadron))])
                   )

    for (obsname, binspace, y_scale) in itertools.product(obsnames,
                                            ['linear', 'log'],
                                            ['linear', 'log']):
        if obsname == 'costheta' and binspace == 'log':
            continue
        # Don't plot lin-log
        if binspace in ['linear', 'lin'] and y_scale == 'log':
            continue

        # Make 2 x 2 figure
        ratio_plot = (obsname in ['costheta', 'z'])
        fig, ax = aesthetic_N_by_M(xlabel=obs_label(obsname),
                             nrows=2, ncols=2,
                             ylabel='EWOC',
                             title=obs_title(obsname),
                             showdate=False,
                             xlim=lims[obsname][(binspace, 'linear')][0],
                             ylim=lims[obsname][(binspace, 'linear')][1],
                             ratio_plot=ratio_plot,
                             rownames=[r"$R_{\rm jet}$"+f"={Rjet}"
                                       for Rjet in kwargs['jet_rad']],
                             colnames=[r"$R_{\rm sub}$"+f"={rsub}"\
                                       for rsub in kwargs['sub_rad']]
                             )

        # Setting up color cycle; parton=green, hadron=gold
        [a.set_prop_cycle((cycler(color=['limegreen', 'gold']))) 
         for a in ax]

        # Iterating over jet and subjet radii:
        # R1 r1, R1 r2, R2 r1, R2 r2
        for i, (Rjet, rsub) in enumerate(
                itertools.product(kwargs['jet_rad'], kwargs['sub_rad'])):
            # Set up this subplot
            for level, hist_dict in zip(['Parton', 'Hadron'],
                        [hist_dicts_parton[i], hist_dicts_hadron[i]]):
                plot_hist_dict(hist_dict, obsname,
                           ax=[ax[2*i], ax[2*i+1]] if ratio_plot
                               else [ax[i]],
                           binspace=binspace,
                           y_scale=y_scale,
                           hist_weight=hist_weight,
                           title=None, label=level)

        # Stamp text and legend
        stamp_ax(ax[0], jet_rad=None, sub_rad=None, obsname=obsname, 
                 **{key: val for key, val in kwargs.items()
                        if key not in ['jet_rad', 'sub_rad']})
        ax[2 if ratio_plot else 1].legend(loc="upper center")

        plt.show()
        

def plot_EWOC_by_frag_temps(load=True, print_every_n=1000,
                            hist_weight=None,
                            **kwargs):
    # Turning the keyword argument for fragmentation temperature into a list
    arg_to_list(kwargs, 'temp')
    # DEBUG
    #input_dict['jet_rad'] = float(input_dict['jet_rad'])
    #input_dict['sub_rad'] = float(input_dict['sub_rad'])

    # Saving data
    for temp in kwargs['temp']:
        save_hist_dict(overwrite=False, print_every_n=print_every_n,
                       **dict(kwargs, temp=temp))

    # Loading and processing data for all given temperatures 
    # (accessed via ```kwargs['temp']```)
    hist_dicts = []
    for temp in kwargs['temp']:
        hist_dicts.append(get_hist_dict(**dict(kwargs, temp=temp),
                                        load=load,
                                        print_every_n=print_every_n))

    # Finding the shared observables in each dict
    obsnames = list(reduce(set.intersection,
                           [set(hist_dict['observables'])
                            for hist_dict in hist_dicts])
                   )

    # Setting up colorbar
    colors, cbar_data = get_colors_colorbar(kwargs['temp'])

    # Plotting
    for obsname in obsnames:
        # yspace -- use itertools
        for binspace, y_scale in zip(['linear', 'log'],
                                     ['linear', 'log']):
            # No logarithmic bins for costheta
            if binspace == 'log' and obsname == 'costheta':
                continue
            # Don't plot lin-log
            if binspace in ['linear', 'lin'] and y_scale == 'log':
                continue


            # Set up figure, axis limits
            fig, ax = aestheticfig(xlabel=obs_label(obsname),
                                   ylabel='EWOC',
                                   title=None,
                                   showdate=False,
                                   xlim=lims[obsname][(binspace, 'linear')][0],
                                   ylim=lims[obsname][(binspace, 'linear')][1],
                                   ratio_plot=(obsname in ['costheta', 'z']))
            # Setting up color cycle
            [a.set_prop_cycle((cycler(color=colors))) for a in ax]

            # Plotting observable EWOC for each subjet radius
            plot_hist_dict(hist_dicts, obsname, ax=ax,
                           binspace=binspace,
                           y_scale=y_scale,
                           hist_weight=hist_weight,
                           **kwargs)

            # Adding subjet radius color bar
            cbar = fig.colorbar(cbar_data, ax=np.array(ax).ravel().tolist())
            cbar.ax.set_ylabel(r'$T_{\rm frag}$')

            # Stamp text
            stamp_ax(ax[0], jet_rad=kwargs['jet_rad'],
                     sub_rad=kwargs['sub_rad'], obsname=obsname, 
                     **{key: val for key, val in kwargs.items()
                            if key not in ['jet_rad', 'sub_rad']})

    plt.show()
