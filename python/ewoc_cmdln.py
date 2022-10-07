import numpy as np
import argparse

# =====================================
# Command Line Parsing 
# =====================================

# Setup
ewoc_parser=argparse.ArgumentParser('Plot EWOC data from Pythia and FastJet')

# ---------------------------------
# Basic Options (Required)
# ---------------------------------
# Event generation
ewoc_parser.add_argument("-n", "--n_events",
            dest="n_events", type=int,
            help="Number of generated events;",
            required=True)
ewoc_parser.add_argument("-l", "--level", 
            dest="qcd_level", type=str,
            help="QCD level of generated events (parton or hadron)",    
            required=True)
ewoc_parser.add_argument("-p", "--process_str", 
            dest="process_str", type=str,
            help="Process for generated events (quark, gluon, qcd, w);",     
            required=True)

# Jet and subjet information
ewoc_parser.add_argument("-j", "--jet_alg", 
            dest="jet_alg_int", type=int,
            help="Jet algorithm (0 [kt], 1 [ca], 2 [antikt]);",
            required=True)
ewoc_parser.add_argument("-s", "--sub_alg", 
            dest="sub_alg_int", type=int,
            help="Subjet algorithm (0 [kt], 1 [ca], 2 [antikt]);",
            required=True)
ewoc_parser.add_argument("-R", "--jet_rad", 
            dest="jet_rad",  # No type, so that it can be a list 
            help="Jet radius (can be 'inf' or 'infty');",
            required=True)
ewoc_parser.add_argument("-r", "--sub_rad", 
            dest="sub_rad",  # No type, so that it can be a list
            help="Subjet radius (or radii);",
            required=True)

# ---------------------------------
# Info for optional variables
# ---------------------------------

default_args    = {'pt_min': 0,
                   'pt_max': np.inf,
                   'E_cm':   4000,
                   'temp': None}

arg_to_str_dict = {'pt_min': 'ptmin',
                   'pt_max': 'ptmax',
                   'E_cm':   'Ecm',
                   'temp':   'temp',}

# ---------------------------------
# Basic Options (Optional)
# ---------------------------------
# Phase Space constraints
ewoc_parser.add_argument("--pt_min", 
            dest="pt_min", type=float,
            help="Minimum value of jet p_T;",
            default=default_args['pt_min'],
            required=False)
ewoc_parser.add_argument("--pt_max", 
            dest="pt_max", type=float,
            help="Maximum value of jet p_T;",
            default=default_args['pt_max'],
            required=False)

ewoc_parser.add_argument("--E_cm", "-E", 
            dest="E_cm", type=float,
            help="Center of mass energy;",
            default=default_args['E_cm'],
            required=False)


# ---------------------------------
# Advanced Options (Optional)
# ---------------------------------
ewoc_parser.add_argument("-T", "--frag_temp", 
            dest="temp",
            help="Temperature of string fragmentation, off by default; see:"\
             +"    * https://pythia.org/latest-manual/Fragmentation.html#anchor25;"\
             +"    * https://arxiv.org/pdf/1610.09818.pdf;",
            default=default_args['temp'],
            required=False)



# =====================================
# Other Command Line Utilities
# =====================================

def arg_to_list(input_dict, key):
    """Turn the argument "key" of an dictionary 
    from a string into a list.
    """
    vals = input_dict[key] 
    try:
        input_dict[key] = [float(val) for val in
                           vals.lstrip('[').rstrip(']').split(',')]

    except AttributeError as e:
        print("The value of the given dictionary at the given key is\
              not a string.")
        print("   ", e)
