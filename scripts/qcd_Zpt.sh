#!/usr/local/bin/fish

# Writing
./write_tools/write_ewocs -n 10000 -l parton -p qcd --write_ewocs false --write_pid_pt 23

# Plotting
./plot_tools/scripts/qcd_Zpt.sh
