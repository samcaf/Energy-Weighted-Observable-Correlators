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
# Basic imports
import os.path
import numpy as np
import matplotlib.pyplot as plt
from math import comb

# Analysis imports
from scipy.signal import find_peaks_cwt

# Plotting imports
from cycler import cycler
from matplotlib import rc
rc('text', usetex=True)

# Iterator utilities
import itertools
from functools import reduce

# Plot utilities
from utils.plot_utils import aestheticfig, aesthetic_N_by_M, stamp
from utils.plot_utils import get_colors_colorbar 
from utils.plot_utils import style_solid, style_dashed, style_scatter, modstyle
from utils.plot_utils import get_hist_edges_centers, text_to_hist

# Physics functions
from utils.qcd_utils import alg_to_string  # jet algorithm names
from utils.qcd_utils import mass_val, scale_name, scale_col
from utils.qcd_utils import eec_nnlo  # eec

# File reading functions
from utils.ewoc_cmdln import ewoc_file_label


# ---------------------------------
# Global Flags
# ---------------------------------
# Make plots with errorbars
plot_errorbars = False
# Make plots with x errorbars (like a histogram)
plot_xerr      = True  # Only applied if plot_errorbars = False 
# Make a scatter plot
plot_scatter   = True  # Only applied if both above are False

# =====================================
# Utilities for reading EWOC data  
# =====================================

# ---------------------------------
# Writing to files 
# ---------------------------------

def save_hist_dict(overwrite, print_every_n=1000,
                  **kwargs):
    file_prefix = ewoc_file_label(**kwargs)
    data_dir = file_prefix + '.txt'
    hist_dir = file_prefix + '_nbins' + str(kwargs['nbins'])\
               + '_normed.npz'

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
    hist_dir = file_prefix + '_nbins' + str(kwargs['nbins'])\
               + '_normed.npz'

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


# =====================================
# Formatting Data 
# =====================================

def ewoc_text_to_dict(filename, pair_obs=None,
                      print_every_n=1000,
                      weighting='energy fraction',
                      run_checks=False):
    """
    Reads in a text file in the form provided by write_ewocs(.cc).
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
        pair_obs = [pair_costheta,  # Usual EEC
                    pair_formtime, pair_mass,  # Helpful observables
                    pair_kt, pair_e1]  # Helpless observables
    if not isinstance(pair_obs, list):
        # Making sure observables are stored in a list
        pair_obs = [pair_obs]
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
        # Additional info for verification:
        num_subjets = -1
        num_subjet_pairs = -1

        # Looping over file
        for i, line in enumerate(file):
            # Read and format the line
            line = line.rstrip("\n")
            info = line.split(" ")

            # - - - - - - - - - - - - - - - - -
            # Misc. information:
            # - - - - - - - - - - - - - - - - -
            if info[0] == "#" or '' in info:
                # Allowing commented or empty lines
                continue

            elif info[0] not in event_headers:
                # Storing additional information from file header
                data_dict[info[0]] = info[-1]
                continue

            # - - - - - - - - - - - - - - - - -
            # Jet Information
            # - - - - - - - - - - - - - - - - -
            if info[0] == "J":
                E_jet = float(info[1])

                if run_checks:
                    # Additional verification that we consider
                    # only subjet pairs within a single jet
                    if num_subjets != -1:
                        # Assert that the number of subjet pairs is
                        # consistent with the number of subjets
                        assert comb(num_subjets, 2) == num_subjet_pairs\
                            or 2.*comb(num_subjets, 2) == num_subjet_pairs,\
                            f"{num_subjets = } not compatible with "+\
                            "number of subjet pairs "+\
                            f"({num_subjet_pairs = }, "+\
                            f"expected {comb(num_subjets, 2) = })."
                    num_subjets = int(info[2])
                    num_subjet_pairs = 0

            # - - - - - - - - - - - - - - - - -
            # Subjet Pair Information
            # - - - - - - - - - - - - - - - - -
            if info[0] == "SP":
                num_subjet_pairs += 1
                # Convert strings to floats
                try:
                    z1, z2, costheta = (float(x) for x in info[1:])
                except ValueError as e:
                    print(e)
                    print(f"    line in file : {info}")
                    print(f"    Next line in file: {file[i+1]}")
                    raise ValueError

                # Store relevant EWOC information

                # Weights
                if weighting == 'energy fraction':
                    data_dict['weights'].append(z1 * z2)
                elif weighting == 'energy total':
                    data_dict['weights'].append(z1 * z2 * E_jet**2)
                else:
                    raise ValueError(f"Invalid weighting {weighting}.")

                # Observables
                for obsname, obs in zip(obsnames, pair_obs):
                    data_dict[obsname].append(obs(E_jet, z1, z2, costheta))

            # - - - - - - - - - - - - - - - - -
            # Event information
            # - - - - - - - - - - - - - - - - -
            if info[0] == "E":
                if i_event%print_every_n == 0:
                    # Better logging
                    print(f"Considered {i_event} events")
                i_event+=1
               
    # After the loop, put weights and observables into numpy arrays
    for obsname in ['weights', *obsnames]:
        data_dict[obsname] = np.array(data_dict[obsname])

    return data_dict


def ewoc_dict_to_hists(data_dict, hist_dict=None,
                       binspaces=['linear', 'log'],
                       **kwargs):
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
        kwargs : dict
            Set of inputs from functions calling this one,
            meant to be received from command line input.
            Used to get number of bins for histograms.

    Returns
    -------
    Dict of histograms containing easy-to-plot ewoc information.
    """
    # Getting number of bins
    print(kwargs)
    nbins = kwargs['nbins']

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
# Pairwise subjet observables (collinear limit) 
# =====================================

# ---------------------------------
# Labels and info
# ---------------------------------

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
    elif observable == 'kt':
        return r'$k_T$ (GeV)'
    elif observable == 'e1':
        return r'$e^{(1)} = k_T/p_T^{(\rm jet)}$'

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
    elif observable == 'kt':
        return r'$k_T$ Subjet EWOC'
    elif observable == 'e1':
        return r'$e^{(1)} = k_T/p_T$ Subjet EWOC'


# ---------------------------------
# Actual observables 
# ---------------------------------

def pair_costheta(Etot, z1, z2, costheta):
    return costheta

def pair_z(Etot, z1, z2, costheta):
    # Not actually used as of 10-24:
    # ctr-f ```if key == 'costheta'``` for more info
    return (1.-costheta)/2.

def pair_m2(Etot, z1, z2, costheta):
    return 2. * Etot**2. * z1 * z2 * (1.-costheta)

def pair_mass(Etot, z1, z2, costheta):
    return np.sqrt(pair_m2(Etot, z1, z2, costheta))

def pair_formtime(Etot, z1, z2, costheta):
    return np.float64(Etot * max(z1, z2)) / pair_m2(Etot, z1, z2, costheta)

def pair_kt(Etot, z1, z2, costheta):
    return Etot * min(z1, z2) * np.sin(np.arccos(costheta))

# Angularities
def pair_angularity(Etot, z1, z2, costheta, beta):
    return min(z1, z2) * np.arccos(costheta)**beta
    # / (z1 + z2)? i.e. angularity of the splitting?

def pair_e1(Etot, z1, z2, costheta):
    return pair_angularity(Etot, z1, z2, costheta, beta=1.0)

# Generalized Jet Energy Correlation Functions
def pair_GECF(Etot, z1, z2, costheta, beta):
    return z1*z2 * np.arccos(costheta)**beta
    # / (z1 + z2)? i.e. GECF of the splitting?

# =====================================
# Plotting Utilites 
# =====================================

# Lower y boundary in log-log plot
_ylog_low = 2e-4
_ylog_high = 2e1

# ---------------------------------
# Dicts and descriptive info 
# ---------------------------------

binedge_dict = {'mass': {'lbin': 1e-5, 'rbin': 2e3},
                'formtime': {'lbin': 5e-4, 'rbin': 1e4}}

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
                      },
         'kt'   : {('linear', 'linear'): 
                        ((0, 1e3), (0, 1.7e-2)),
                      ('linear', 'log'): 
                        ((0, 150), (1e-4, .12)),
                      ('log', 'linear'):
                        ((2e-3, 1e3), (0, .9)),
                      ('log', 'log'):
                        ((1.0, 1e3), (1e-4, .7))
                      },
         'e1'       : {('linear', 'linear'): 
                        ((0, .25), (0, 20)),
                      ('linear', 'log'): 
                        ((0, .25), (.1, 50)),
                      ('log', 'linear'):
                        ((3e-5, 1e0), (0, .5)),
                      ('log', 'log'):
                        ((1e-6, .2), (.1, 50))
                      },
        }


def expected_obs_peak(obsname, **kwargs):
    E_cm = kwargs['E_cm']
    mass = mass_val[kwargs['process_str']]

    if obsname == 'mass':
        peak = mass
    elif obsname == 'costheta':
        peak = 1 - 8.*mass**2./E_cm**2.
    elif obsname == 'z':
        peak = (2*mass/E_cm)**2.
    elif obsname == 'formtime':
        expected_E_hard = E_cm/3
        peak = expected_E_hard/mass**2.
    elif obsname == 'kt':
        expected_E_soft = E_cm/10
        peak = expected_E_soft * 2.*mass/E_cm
    elif obsname == 'e1':
        expected_E_soft = E_cm/10
        peak = expected_E_soft * 2.*mass/E_cm**2
    else:
        raise AssertionError(f"Invalid {obsname = }")

    return peak


def expected_peaks_scales_colors(obsname, **kwargs):
    exp_peaks = [expected_obs_peak(obsname, **kwargs)]
    scales = [scale_name[kwargs['process_str']]]
    colors = [scale_col[kwargs['process_str']]]

    # If we expect additional peaks from other physics
    if "Z" in kwargs['s_channel']:
        exp_peaks.append(expected_obs_peak(obsname,
                            **dict(kwargs, process_str='z')))
        scales.append(scale_name['z'])
        colors.append(scale_col['z'])
    if kwargs['process_str'] == 'top':
        exp_peaks.append(expected_obs_peak(obsname,
                            **dict(kwargs, process_str='w')))
        scales.append(scale_name['w'])
        colors.append(scale_col['w'])
    
    return exp_peaks, scales, colors


def hist_title(observable, jet_radius, subjet_alg):
    """Title for a histogram plotting the EWOC for the given
    subjet pair property, jet radius, and subjet algorithm.
    """
    # radius info
    title = r'EWOC, $R_{\rm jet}$=' + f'{float(jet_radius):.1f}'

    # Add subjet alg info
    title += f', {alg_to_string(subjet_alg)} Subjets'

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
                  f"({alg_to_string(kwargs['jet_alg'],latex=False).upper()}"\
                  +(f" {jet_rad}" if jet_rad is not None else '') + ", "\
                  +f"{alg_to_string(kwargs['sub_alg'],latex=False).lower()}"\
                  +(f" {sub_rad}" if sub_rad is not None else '')\
                  +fr") \textbf{{{obs_title(obsname)}}}",
               'line_1':
                  f"{kwargs['n_events']} {kwargs['E_cm']/1000} TeV "
                  +(f"{kwargs['process_str'].upper()} events, "\
                  if kwargs['process_str'] != 'top' else "Top events, ")
                  +r"\texttt{Pythia 8.307}",
               'line_2':
                  str(kwargs['pt_min']) + r" $< p_T <$ "\
                  + str(kwargs['pt_max'])
            }

    # Make the stamp!
    stamp(*stamp_loc, ax, **stamp_text)


# ---------------------------------
# General plotting utils 
# ---------------------------------

def plot_hist_dict(hist_dicts, obsname,
                   ax=None,
                   binspace=None, y_scale='linear', 
                   xlim='default',
                   ylim='default',
                   title=None,
                   show_analytic=True,
                   hist_weight=None,
                   label=None,
                   plot_expectation=True,
                   **kwargs):
    # Showing analytic EWOC for only specific observables
    show_analytic = obsname in ['costheta', 'z'] and show_analytic

    # Extra weighting factors for plotting the given histogram
    if hist_weight is None:
        # Default: no additional weight
        hist_weight = lambda x: np.ones(len(x))

    # Figure setup
    if ax is None:
        _, ax = aestheticfig(xlabel=obs_label(obsname),
                             ylabel='EWOC',
                             title=None,
                             showdate=False,
                             xlim=xlim if xlim != 'default' else
                                  lims[obsname][(binspace, 'linear')][0],
                             ylim=ylim if ylim != 'default' else
                                  lims[obsname][(binspace, 'linear')][1],
                             ratio_plot=(obsname in ['costheta', 'z']))

    # Plot line showing an order-of-magnitude expectation for EWOC peaks
    if plot_expectation:
        exp_peaks, scales, colors = expected_peaks_scales_colors(obsname,
                                                                 **kwargs)

        # Plotting vertical lines
        ax[0].vlines(exp_peaks, -1, 1 if y_scale == 'log' else 1e7, colors=colors)
        if show_analytic:
            ax[1].vlines(exp_peaks, -1, 1e7, colors=colors)
        # Adding text 
        ylow, yhigh = ax[0].get_ylim()
        ymid = (ylow+yhigh)/2. if y_scale in ['lin', 'linear']\
            else np.sqrt((ylow+2*_ylog_low)*yhigh)
        for iscale, scale in enumerate(scales):
            scale_frac = (iscale+1/2.)/len(scales)
            if y_scale in ['lin', 'linear']:
                y_text = (1.-scale_frac)*ylow + scale_frac*yhigh
            else:
                y_text = (2*_ylog_low+ylow)**(1.-scale_frac) * ymid**scale_frac
            ax[0].text(exp_peaks[iscale], y_text,
                       ("Set by " if iscale==0 else '') + scale,
                       rotation=90, verticalalignment='center')

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
        elif plot_xerr:
            ax[0].errorbar(xs, ys, xerr=(xs - x_edges[:-1], x_edges[1:] - xs),
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
            elif plot_xerr:
                ax[1].errorbar(xs, ys/analytic,
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
        ax[0].set_ylim(_ylog_low, _ylog_high)
        ax[0].set_yscale('log')

    return


def plot_hist_dict_peaks(hist_dicts, obsname,
                         x_var=None, xlabel=None,
                         binspace='log',
                         x_scale='lin', y_scale='log', 
                         title=None,
                         plot_expectation=True,
                         **kwargs):
    # x-axis variable setup (default r/R)
    if x_var is None:
        def x_var(sub_rad, jet_rad, **kwargs):
            return float(sub_rad)/float(jet_rad)
        xlabel = r"$r_{\rm sub} / R_{\rm jet}$"

    # Figure setup
    fig, ax = aestheticfig(xlabel=xlabel,
                ylabel=obs_title(obsname) + " Peak")

    # Plotting for each hist dict
    for hist_dict in hist_dicts:
        # Getting the values on the x axis for this dict
        x_val = x_var(**hist_dict)
        hist, _, bin_centers = hist_dict[obsname][binspace]
        # for peak_ind in find_peaks(hist)[0]:
        for peak_ind in find_peaks_cwt(hist, np.arange(1,10)):
            ax[0].scatter(x_val, bin_centers[peak_ind],
                          color='cornflowerblue')

    if x_scale == 'log':
        [a.set_xscale('log') for a in ax]
    if y_scale == 'log':
        ax[0].set_ylim(_ylog_low, _ylog_high)
        ax[0].set_yscale('log')

    # Plotting what we might expect for the peak value
    if plot_expectation:
        exp_peaks, scales, colors = expected_peaks_scales_colors(obsname,
                                                                 **kwargs)
        # Plotting vertical lines
        ax[0].hlines(exp_peaks, -1, 1e7, colors=colors)
        # Adding text 
        xlow, xhigh = ax[0].get_xlim()
        xmid = (xlow+xhigh)/2. if x_scale in ['lin', 'linear']\
            else np.sqrt((xlow+2e-8)*xhigh)
        for iscale, scale in enumerate(scales):
            scale_frac = (iscale+1/2.)/len(scales)
            if x_scale in ['lin', 'linear']:
                x_text = (1.-scale_frac)*xlow + scale_frac*xhigh
            else:
                x_text = (2e-8+xlow)**(1.-scale_frac) * xmid**scale_frac
            ax[0].text(exp_peaks[iscale], x_text,
                       ("Set by " if iscale==0 else '') + scale,
                       rotation=0, verticalalignment='center')



# ---------------------------------
# Plotting by subjet radius 
# ---------------------------------

def plot_EWOC_by_sub_rads(load=True, print_every_n=1000,
                          hist_weight=None,
                          obsnames=None,
                          overwrite=False,
                          show=True,
                          **kwargs):
    # Typesetting keyword arguments
    for key, value in kwargs.items():
        if key == 'sub_rad':
            if not isinstance(value, list):
                kwargs[key] = [value]
        else:
            if isinstance(value, list):
                assert len(value) == 1, "Received too many values"\
                    +f"({len(value)} for {key}."
                kwargs[key] = value[0]
            
    if not isinstance(kwargs['sub_rad'], list):
        kwargs['sub_rad'] = [kwargs['sub_rad']]
            
    # Saving data
    for rsub in kwargs['sub_rad']:
        save_hist_dict(overwrite=overwrite, print_every_n=print_every_n,
                       **dict(kwargs, sub_rad=rsub))
    
    # Loading and processing data all given subjet radii
    # (accessed via ```kwargs['sub_rad']```)
    hist_dicts = []

    for rsub in kwargs['sub_rad']:
        hist_dicts.append(get_hist_dict(**dict(kwargs, sub_rad=rsub),
                                        load=load,
                                        print_every_n=print_every_n))

    # Finding the shared observables in each dict
    if obsnames is None:
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

            # Adding subjet radius color bar
            cbar = fig.colorbar(cbar_data, ax=np.array(ax).ravel().tolist())
            cbar.ax.set_ylabel(r'$R_{\rm sub}$')

            # Stamp text
            stamp_ax(ax[0], jet_rad=kwargs['jet_rad'],
                     sub_rad=None, obsname=obsname,
                     **{key: val for key, val in kwargs.items()
                            if key not in ['jet_rad', 'sub_rad']})

        # Plotting the peaks of these distributions
        if obsname != 'costheta':
            plot_hist_dict_peaks(hist_dicts, obsname, binspace='log', **kwargs)

    if show:
        plt.show()

    return


# ---------------------------------
# Plotting parton vs. hadron by (sub)jet radius 
# ---------------------------------

def plot_EWOC_pvh_by_rads(load=True, print_every_n=1000,
                          hist_weight=None,
                          obsnames=None,
                          overwrite=False,
                          show=True,
                          **kwargs):
    # Turning the keyword argument for (sub)jet radii into a list
    for key, value in kwargs.items():
        if key in ['sub_rad', 'jet_rad']:
            assert len(value) == 2, "Expected two radii"\
                + f"of type {rad_type} (received {len(value)})."
        else:
            if isinstance(value, list):
                assert len(value) == 1, "Received too many values"\
                    +f"({len(value)} for {key}."
                kwargs[key] = value[0]

    # Saving data
    for (Rjet, rsub, level) in itertools.product(
                                kwargs['jet_rad'],
                                kwargs['sub_rad'],
                                ['parton', 'hadron']):
        save_hist_dict(overwrite=overwrite, print_every_n=print_every_n,
                       **dict(kwargs, qcd_level=level,
                              jet_rad=Rjet, sub_rad=rsub))

    # Getting observable dictionaries for all given subjet radii
    # (accessed via ```kwargs['sub_rad']```)
    hist_dicts_parton = []
    hist_dicts_hadron = []
      
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
    if obsnames is None:
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

        # Setting up color cycle; parton=blue, hadron=red/brown
        [a.set_prop_cycle((cycler(color=['royalblue', 'chocolate']))) 
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
                           title=None, label=level,
                           plot_expectation=(level=='Parton'), 
                           **kwargs)

        # Stamp text and legend
        stamp_ax(ax[0], jet_rad=None, sub_rad=None, obsname=obsname, 
                 **{key: val for key, val in kwargs.items()
                        if key not in ['jet_rad', 'sub_rad']})
        ax[2 if ratio_plot else 1].legend(loc="upper center")

        if show:
            plt.show()
        

# ---------------------------------
# Plotting by fragmentation temperature
# ---------------------------------

def plot_EWOC_by_frag_temps(load=True, print_every_n=1000,
                            hist_weight=None,
                            obsnames=None,
                            overwrite=False,
                            show=True,
                            **kwargs):
    # Turning the keyword argument for fragmentation temperature into a list
    for key, value in kwargs.items():
        if key == 'temp':
             if not isinstance(value, list):
                kwargs[key] = [value]
        else:
            if isinstance(value, list):
                assert len(value) == 1, "Received too many values"\
                    +f"({len(value)} for {key}."
                kwargs[key] = value[0]
            
    # Saving data
    for temp in kwargs['temp']:
        save_hist_dict(overwrite=overwrite, print_every_n=print_every_n,
                       **dict(kwargs, temp=temp))

    # Loading and processing data for all given temperatures 
    # (accessed via ```kwargs['temp']```)
    hist_dicts = []
    for temp in kwargs['temp']:
        hist_dicts.append(get_hist_dict(**dict(kwargs, temp=temp),
                                        load=load,
                                        print_every_n=print_every_n))

    # Finding the shared observables in each dict
    if obsnames is None:
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

    if show:
        plt.show()
