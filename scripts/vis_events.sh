#!/usr/bin/env bash
# Script for quick visualization
# Example usage:
# ./scripts/vis_events.sh -l hadron -p qcd -j akt -s ca -R 2 -r 1

# Writing an event
./write_tools/write_ewocs --write_event true --write_ewocs false $@ -n 1
# Note: default looks at the first event generated by pythia, but adding the option
#       ```-n [N]```
#       will visualize the Nth event

# write_ewocs stores a pointer to the filename with event info in another file...
value=`cat event_vis_pointer.txt`
echo "Visualizing event stored at "$value

# Visualizing
./plot_tools/event_vis --filename $value scatter_vis

echo "Erasing the evidence >:)"
rm $value
rm event_vis_pointer.txt
