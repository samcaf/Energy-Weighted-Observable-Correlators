#!/usr/local/bin/fish

# Default arguments
set -Ux generate_events false
set -Ux overwrite_hists false
set -Ux show_plots true

# Checking for user input
getopts $argv | while read -l key value
    switch $key
        case generate_events
            set -Ux generate_events $value
        case overwrite_hists
            set -Ux overwrite_hists $value
        case show_plots
            set -Ux show_plots $value
    end
end
