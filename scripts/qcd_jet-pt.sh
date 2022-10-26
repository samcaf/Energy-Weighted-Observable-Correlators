#!/usr/local/bin/fish

# Writing
./write_tools/write_ewocs -n 10000 -l parton -p qcd --pt_min 50 --pt_max 3000 --write_ewocs false --write_jet_pt true

# Plotting
./plot_tools/plot_ewocs -n 10000 -l parton -p qcd --pt_min 50 --pt_max 3000 --plot_jet_pt
