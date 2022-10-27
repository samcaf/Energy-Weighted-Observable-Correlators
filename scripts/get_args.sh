#!/usr/local/bin/fish

# =======================
# Default arguments
# =======================
# Use e+e- algorithms when finding jets
set -Ux ee_algs true
# Generate events with an event generator
set -Ux generate_events false
# Overwrite existing histogram data
# (i.e. re-analyze event generator data)
set -Ux overwrite_hists false
# Show plots
set -Ux show_plots true

# =======================
# Checking for user input
# =======================
getopts $argv | while read -l key value
    switch $key
        case ee_algs
            set -Ux ee_algs $value
        case generate_events
            set -Ux generate_events $value
        case overwrite_hists
            set -Ux overwrite_hists $value
        case show_plots
            set -Ux show_plots $value
    end
end

# If we are told to use e+e- algorithms, change the scripts (changed back by `./clear_args.sh`)
if $ee_algs
    echo "Setting up e+ e- algorithms in scripts:"
    find scripts/ -type f -print0 | xargs -0 sed -i '' -E 's|' kt'|' ee_kt'|g'
    find scripts/ -type f -print0 | xargs -0 sed -i '' -E 's|' akt'|' ee_akt'|g'
    find scripts/ -type f -print0 | xargs -0 sed -i '' -E 's|' ca'|' ee_ca'|g'
end
