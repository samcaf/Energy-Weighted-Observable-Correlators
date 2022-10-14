import numpy as np
import matplotlib.pyplot as plt


def ptyphis_from_jet(ptyphis, ijet):
    return np.array([ptyphi for i, ptyphi
                    in enumerate(ptyphis[:,1:])
                    if event[i][0] == ijet])


def ptyphis_with_user_index(ptyphis, user_index):
    return np.array([ptyphi for i, ptyphi
                    in enumerate(ptyphis[:,:-1])
                    if event[i][-1] == user_index])


def event_vis_file(): -> str
    return "output/qcd_parton/jetR1-0/aktjet/subR0-00_casub_10000evts_ptmin50-0_ptmax3000-0_Ecm10000-0.txt_vis-evt0.txt"


def visualize_ptyphis(ptyphis, ax, color):
    pts, ys, phis = ptyphis[:,0], ptyphis[:,1], ptyphis[:,2]

    hist, y_edges, phi_edges = np.histogram2d(ys, phis, bins=(20,20), weights=pts)
    y_pos, phi_pos = np.meshgrid(y_edges[:-1]+y_edges[1:], phi_edges[:-1]+phi_edges[1:])

    y_pos = y_pos.flatten()/2.
    phi_pos = phi_pos.flatten()/2.
    hist_pos = np.zeros_like (y_pos)

    dy = y_edges[1] - y_edges[0]
    dphi = phi_edges[1] - phi_edges[0]
    dz = hist.flatten()

    ax.bar3d(y_pos, phi_pos, hist_pos, dy, dphi, dz, color=color, zsort='average')


def visualize_event(event_ptyphis, colors=None):
    if colors is None:
        colors = ['orangered', 'forestgreen', 'cornflowerblue']

    # Finding ptyphis associated with particles and ghosts
    part_ptyphis = ptyphis_with_user_index(event_ptyphis, 'P')
    ghost_ptyphis = ptyphis_with_user_index(event_ptyphis, 'G')

    # Setting up figure with a 3D canvas
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Looping over all jets
    ijet = 0
    while True:
        ijet += 1
        # All pseudojets associated with this jet
        this_jet_particles = ptyphis_from_jet(part_ptyphis, ijet)
        this_jet_ghosts = ptyphis_from_jet(ghost_ptyphis, ijet)

        visualize_ptyphis(this_jet_particles, colors[ijet])
        visualize_ptyphis(this_jet_ghosts, colors[ijet])
        # scale_lightness(colors[ijet], .3)

    # Finishing figure details
    plt.title("Event Visualization")
    plt.xlabel(r"Rapidity ($y$)")
    plt.ylabel(r"Azimuthal Angle ($\phi$)")
    plt.savefig("this_evt_vis.pdf")
    plt.show()

