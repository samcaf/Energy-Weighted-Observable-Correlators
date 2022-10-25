#!/usr/local/bin/fish

# ==========================
# Getting arguments
# ==========================
./scripts/get_args.sh $argv


# ==========================
# Writing
# ==========================
if $generate_events
for rsub in 0.0 0.1 0.2 0.3 0.4 0.5
  ./write_tools/write_ewocs -n 10000 -l parton -p qcd -j akt -s ca -R 1.0 -r $rsub --pt_min 50 --pt_max 3000
end
end


# ==========================
# Plotting
# ==========================
if begin ; $overwrite_hists ; and $show_plots ; end
  ./plot_tools/plot_ewocs -n 10000 -l parton -p qcd -j akt -s ca -R 1.0 -r 0 0.1 .2 .3 .4 .5 --pt_min 50 --pt_max 3000 --plot_type sub_rads --overwrite_hists true --show_plots true
else if $overwrite_hists
  ./plot_tools/plot_ewocs -n 10000 -l parton -p qcd -j akt -s ca -R 1.0 -r 0 0.1 .2 .3 .4 .5 --pt_min 50 --pt_max 3000 --plot_type sub_rads --overwrite_hists true --show_plots false
else if $show_plots
  ./plot_tools/plot_ewocs -n 10000 -l parton -p qcd -j akt -s ca -R 1.0 -r 0 0.1 .2 .3 .4 .5 --pt_min 50 --pt_max 3000 --plot_type sub_rads --overwrite_hists false --show_plots true
end


# ==========================
# Clearing Arguments
# ==========================
set --erase generate_events
set --erase overwrite_hists
set --erase show_plots
