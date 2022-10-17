#!/usr/bin/env bash
# Script for quick visualization

# Writing an event
./write_tools/write_ewocs -n 1 --write_event true --write_ewocs false $@

# write_ewocs stores a pointer to the filename with event info in another file...
value=`cat event_vis_pointer.txt`
echo "Visualizing event stored at "$value;

# Visualizing
./plot_tools/event_vis --filename $value scatter_vis

# Erasing the evidence
rm $value
rm event_vis_pointer.txt
