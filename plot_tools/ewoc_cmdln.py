import argparse

from numpy import inf
from matplotlib.pyplot import show

from sys import exit

from plot_utils import aestheticfig, text_to_hist

from qcd_utils import alg_to_string


# =====================================
# EWOC filenames 
# =====================================

def process_folder(qcd_level, process_str, **kwargs):
    """
    Returns a folder used to store process information given
    a set of command line inputs

    Parameters
    ----------
        Command line arguments; see ewoc_cmdln.py

    Returns
    -------
        string : The relevant folder.

    """
    proc_folder = "output/" + process_str + "_" + qcd_level
    if kwargs['s_channel'] != default_args['s_channel']:
        proc_folder += "_schan_" + kwargs['s_channel']
    proc_folder += "/"
    return proc_folder


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
    return process_folder(qcd_level, process_str, **kwargs)\
        + f"jetR{float(jet_rad):.1f}/".replace(".", "-")\
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


# =====================================
# Command Line Parsing 
# =====================================

# Setup
ewoc_parser=argparse.ArgumentParser('Plot EWOC data from Pythia and FastJet')

# ---------------------------------
# Additional Actions
# ---------------------------------

# Action class for reading jet algorithm information
class parse_jet_alg(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values in ['akt', 'antikt', '2', 2]:
            alg_int = 2
        elif values in ['ca', '1', 1]:
            alg_int = 1
        elif values in ['kt', '0', 0]:
            alg_int = 0
        else:
            raise AttributeError("Invalid jet algorithm value: {values}")
        setattr(namespace, self.dest, alg_int)


# Accept either single arguments or list of arguments
class store_val_or_list(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        val_or_list = [float(val) for val in
                       values.lstrip('[').rstrip(']').split(',')]
        if len(val_or_list) == 1:
            val_or_list = val_or_list[0]
        setattr(namespace, self.dest, val_or_list)


# Plot pT associated with a particular PID
class plot_pid_pT(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        pt_filename = process_folder(**vars(namespace))\
            + str(values) + "-pt-spectrum.txt" 
        text_to_hist(pt_filename, xlabel=fr"$p^{{\rm PID = {values}}}_T$")
        show()
        exit()


# Plot fraction of energy given to softer subjet
class plot_soft_energy(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        kwargs = vars(namespace) 

        for rad_type in ['sub_rad', 'jet_rad']:
            if not isinstance(kwargs[rad_type], list):
                kwargs[rad_type] = [kwargs[rad_type]]

        fig, ax = aestheticfig(ratio_plot=False, showdate=False,
                               xlabel=r"$E_{\rm soft\ subjet}$")

        def soft_energy(z1, z2, Etot):
            return Etot * min(z1, z2)

        for jet_rad in kwargs['jet_rad']:
            for sub_rad in kwargs['sub_rad']:
                filename = ewoc_file_label(**dict(kwargs,
                                          jet_rad=jet_rad,
                                          sub_rad=sub_rad))
                text_to_hist(filename, use_cols=[1, 2, 3],
                             col_func=soft_energy,
                             fig=fig, ax=ax)
        show()
        exit()


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
            action=parse_jet_alg, dest="jet_alg_int",
            default=2,
            help="Jet algorithm: 0 [kt], 1 [ca], or "\
                 "2 [akt, antikt] (default);")
ewoc_parser.add_argument("-s", "--sub_alg", 
            action=parse_jet_alg, dest="sub_alg_int",
            default=1,
            help="Subjet algorithm (0 [kt], 1 [ca] (default), "\
                 +"or 2 [akt, antikt];")

ewoc_parser.add_argument("-R", "--jet_rad", 
            action=store_val_or_list,
            dest="jet_rad", 
            help="Jet radius;", default=1.0)
ewoc_parser.add_argument("-r", "--sub_rad", 
            action=store_val_or_list, dest="sub_rad",
            help="Subjet radius (or radii);", default=0.1)

# ---------------------------------
# Info for optional variables
# ---------------------------------

default_args    = {'pt_min':    0,
                   'pt_max':    inf,
                   'E_cm':      4000,
                   'temp':      None,
                   's_channel': 'gmZ'}

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

ewoc_parser.add_argument("--s_channel", 
            dest="s_channel", type=str,
            help="Allowed s-channel mediators;",
            default=default_args['s_channel'],
            required=False)


# ---------------------------------
# Advanced Options (Optional)
# ---------------------------------
ewoc_parser.add_argument("-T", "--frag_temp", 
            dest="temp",
            action=store_val_or_list,
            help="Temperature of string fragmentation, off by default; see:"\
             +"\n    * https://pythia.org/latest-manual/Fragmentation.html#anchor25;"\
             +"\n    * https://arxiv.org/pdf/1610.09818.pdf;",
            default=default_args['temp'],
            required=False)


# ---------------------------------
# Easy Plotting Options (Optional)
# ---------------------------------

ewoc_parser.add_argument("--plot_type",
            dest="plot_type", type=str,
            help="Argument which sets the type of plot to create."\
                 +"Can be:\n*'sub_rad';\n'pvh';\n*'thermal'.",
            default='sub_rad',
            required=False)
                
ewoc_parser.add_argument("--plot_pt_pid", 
            action=plot_pid_pT,
            help="Optional argument which histograms the pT of particles "
                 +"with given particle ID that appear in the given "
                 +"process."\
                 +"\nRelies on input from the"\
                 +"```--write_pt_pid``` option of "
                 +"```write_ewocs```.",
            default=None,
            required=False)

ewoc_parser.add_argument("--plot_soft_energy", 
            action=plot_soft_energy,
            help="Optional argument which histograms the energy "\
                  "carried by the softer subjet in every pair of"\
                  "subjets found using the specified jet parameters.",
            default=None,
            required=False)

