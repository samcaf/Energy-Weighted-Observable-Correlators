import argparse

from matplotlib.pyplot import legend, show

from sys import exit

from utils.plot_utils import aestheticfig, text_to_hist

from utils.qcd_utils import alg_to_string


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
                         + f"{float(kwargs[key]):.1f}".replace(".", "-")

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


# Plot pT associated with a particular PID
class plot_pid_pT(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        pt_filename = process_folder(**vars(namespace))\
            + str(values) + "-pt-spectrum.txt" 
        text_to_hist(pt_filename,
                     xlabel=fr"$p^{{\rm PID = {values}}}_T$",
                     nbins=2500, binspace='lin')
                     #nbins=25, binspace='log')
        if vars(namespace)['show']:
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
            z1, z2, Etot = float(z1), float(z2), float(Etot)
            return Etot * min(z1, z2)

        def valid_row(row_info):
            return row_info[0] == 'SP'

        for jet_rad in kwargs['jet_rad']:
            for sub_rad in kwargs['sub_rad']:
                filename = ewoc_file_label(**dict(kwargs,
                                          jet_rad=jet_rad,
                                          sub_rad=sub_rad))
                filename += '.txt'
                text_to_hist(filename, use_cols=[1, 2, 3],
                             col_func=soft_energy,
                             use_rows=valid_row,
                             fig=fig, ax=ax)
        if vars(namespace)['show']:
            show()
        exit()


# Plot pT associated with jets 
class plot_jet_pT(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        fig, ax = aestheticfig(ratio_plot=False, showdate=False,
                               xlabel="Jet Property (GeV)")

        # Finding pT spectrum directly
        pt_filename = ewoc_folder(**vars(namespace))\
            + "jet-pt-spectrum"\
            + "_min"+str(round(vars(namespace)['pt_min']))\
            + "_max"+str(round(vars(namespace)['pt_max']))\
            + ".txt" 
        text_to_hist(pt_filename,
                     xlabel=r"$p^{\rm jet}_T$",
                     nbins=25, binspace='lin',
                     upper_x=2200,
                     label=r"$p_T$ Spectrum",
                     color="steelblue",
                     fig=fig, ax=ax)

        # As a common sense check, finding it independently
        # jet info in the relevant EWOC file
        def valid_row(row_info):
            return row_info[0] == 'J'

        filename = ewoc_file_label(**vars(namespace))
        filename += '.txt'

        text_to_hist(filename, use_cols=[1],
                     use_rows=valid_row,
                     nbins=25, binspace='lin',
                     upper_x=2200,
                     label=r"Energy Spectrum",
                     color="lightcoral",
                     fig=fig, ax=ax)
        legend()

        if vars(namespace)['show']:
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
            nargs='+', type=float,
            dest="jet_rad", 
            help="Jet radius;", default=1.0)
ewoc_parser.add_argument("-r", "--sub_rad", 
            nargs='+', type=float,
            dest="sub_rad",
            help="Subjet radius (or radii);", default=0.1)

# ---------------------------------
# Info for optional variables
# ---------------------------------

default_args    = {'pt_min':    0,
                   'pt_max':    1e5,
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
            dest="temp", type=float, nargs='+',
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
                 +"Can be:\n*'sub_rad';\n'pvh';\n*'thermal';",
            default='sub_rad',
            required=False)
                
ewoc_parser.add_argument("--plot_pid_pt", 
            action=plot_pid_pT,
            help="Optional argument which histograms the pT of particles "
                 +"with given particle ID that appear in the given "
                 +"process."\
                 +"\nRelies on input from the"\
                 +"```--write_pid_pt``` option of "
                 +"```write_ewocs```;",
            default=None,
            required=False)

ewoc_parser.add_argument("--plot_jet_pt", nargs=0, 
            action=plot_jet_pT,
            help="Optional argument which histograms the pT of jets "
                 +"with given parameters in the given process."\
                 +"\nRelies on input from the"\
                 +"```--write_jet_pt``` option of "
                 +"```write_ewocs```;",
            required=False)

ewoc_parser.add_argument("--plot_soft_energy", 
            action=plot_soft_energy, nargs=0,
            help="Optional argument which histograms the energy "\
                  "carried by the softer subjet in every pair of"\
                  "subjets found using the specified jet parameters;",
            default=None,
            required=False)


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

ewoc_parser.add_argument("--overwrite_hists",
            type=str2bool, nargs='?',
            const=True, default=False,
            dest="overwrite", 
            help="Flag determining whether to overwrite existing data"
                +" when generating plots;")
ewoc_parser.set_defaults(feature=True)

ewoc_parser.add_argument("--show_plots", 
            type=str2bool, nargs='?',
            const=True, default=True,
            dest="show", 
            help="Flag determining whether to show plots;")
