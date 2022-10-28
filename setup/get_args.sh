#!/usr/local/bin/fish

# =======================
# Default arguments
# =======================
# Use e+e- algorithms when finding jets
set -Ux alg_prefix ""
# Generate events with an event generator
set -Ux generate_events false
# Overwrite existing histogram data
# (i.e. re-analyze event generator data)
set -Ux overwrite_hists false
# with some number of bins (default 100)
set -Ux nbins 100
# Show plots
set -Ux show_plots true

# =======================
# Checking for user input
# =======================
getopts $argv | while read -l key value
    switch $key
        # Option for ee jet algorithms
        case ee ee_algs ee_algorithms
            set -Ux alg_prefix "ee_"
        # Option to generate new events
        case generate_events
            set -Ux generate_events $value
        # Option to generate new histograms
        case overwrite_hists
            set -Ux overwrite_hists $value
        # Option for number of histogram bins
        case nbins n_bins
            set -Ux nbins $value
        # Option to show plots
        case show_plots
            set -Ux show_plots $value
    end
end
