#!/usr/local/bin/fish

./setup/get_args.sh $argv


# Writing
if $generate_events
  ./write_tools/write_ewocs -n 10000 -l parton -p qcd --pt_min 50 --pt_max 3000 --write_ewocs true --write_jet_pt true -v 2
end

# Plotting
./plot_tools/plot_ewocs --nbins $nbins -n 10000 -l parton -p qcd --pt_min 50 --pt_max 3000 --plot_jet_pt

./setup/clear_args.sh
