import os.path
import numpy as np

from ewoc_cmdln import ewoc_file_label

# =====================================
# Utilities for reading EWOC data  
# =====================================

# ---------------------------------
# Writing to files 
# ---------------------------------

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


# =====================================
# Formatting Data 
# =====================================

def ewoc_text_to_dict(filename, pair_obs=None,
                      print_every_n=1000):
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


