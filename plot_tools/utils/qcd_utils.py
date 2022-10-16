import numpy as np

# =====================================
# QCD and Analytic Expectations
# =====================================

# ---------------------------------
# Mass scales 
# ---------------------------------

# Mass Scales (in GeV)
M_Z = 91.19
M_W = 80.38
M_t = 172.76
Lambda_QCD = .245

mass_val = {'qcd' : Lambda_QCD,
             'z'   : M_Z,
             'w'   : M_W,
             'top' : M_t}

scale_name = {'qcd' : r"$\Lambda_{\rm QCD}$",
              'z'   : r"$m_Z$",
              'w'   : r"$m_W$",
              'top' : r"$m_t$"}

scale_col = {'qcd' : 'rebeccapurple',
             'z'   : 'cadetblue',
             'w'   : 'plum',
             'top' : 'peru'}


# ---------------------------------
# QCD 
# ---------------------------------

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



# ---------------------------------
# EECs 
# ---------------------------------

def eec_nnlo(z, mu):
    """NNLO EEC."""
    a_s = alpha_s(mu) / (4 * np.pi)
    # Not quite right yet, need to run at higher acccuracy

    z = z.astype(float)

    nlo_piece  = -11.5333*np.log(z) + 81.4809
    nnlo_piece = 45.1489*np.log(z)**2. - 1037.73*np.log(z) + 2871.36

    return (2.*a_s + a_s**2.*nlo_piece + a_s**3.*nnlo_piece)/(z*(1-z))


def plot_eec_analytic(ax, observable, binspace, Q=2000, acc='nnlo'):
    """Plots an analytic EEC to the given accuracy on the given axes."""
    assert acc == 'nnlo', "Invalid accuracy " + str(acc)

    # Plot the full analytic EEC
    if observable in ['z', 'zs']:
        xs = np.concatenate((np.linspace(1e-8, 1, 250),
                             np.logspace(-8, 0, 250)))
        xs = np.sort(xs)
        ax.plot(xs, eec_nnlo(xs, Q), 'k--',
                linewidth=2, zorder=4)
    if observable in ['cos', 'costhetas']:
        xs = np.linspace(-1, 1, 500)
        ax.plot(xs, eec_nnlo((1-xs)/2, Q), 'k--',
                linewidth=2, zorder=4)
    return


# ---------------------------------
# Jet algorithm utils
# ---------------------------------

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
