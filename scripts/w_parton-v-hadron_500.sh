#!/usr/local/bin/fish

# ==========================
# Getting arguments
# ==========================
./scripts/get_args.sh $argv


# ==========================
# Writing
# ==========================
if $generate_events
for Rjet in 1000.0 1.0;
  for rsub in 0.0 0.1;
    # kt with -s 1, ca with -s 0
    ./write_tools/write_ewocs -n 10000 -l parton -p w -j akt -s ca -R $Rjet -r $rsub --pt_min 50 --pt_max 3000 -E 500;
    ./write_tools/write_ewocs -n 10000 -l hadron -p w -j akt -s ca -R $Rjet -r $rsub --pt_min 50 --pt_max 3000 -E 500;
  end;
end;
end

# ==========================
# Plotting
# ==========================
if begin ; $overwrite_hists ; and $show_plots ; end
  ./plot_tools/plot_ewocs -n 10000 -l parton -p w -j akt -s ca -R 1000 1 -r 0 .1 --pt_min 50 --pt_max 3000 -E 500 --plot_type pvh --overwrite_hists true --show_plots true
else if $overwrite_hists
  ./plot_tools/plot_ewocs -n 10000 -l parton -p w -j akt -s ca -R 1000 1 -r 0 .1 --pt_min 50 --pt_max 3000 -E 500 --plot_type pvh --overwrite_hists true --show_plots false
else if $show_plots
  ./plot_tools/plot_ewocs -n 10000 -l parton -p w -j akt -s ca -R 1000 1 -r 0 .1 --pt_min 50 --pt_max 3000 -E 500 --plot_type pvh --overwrite_hists false --show_plots true
end


# ==========================
# Clearing Arguments
# ==========================
set --erase generate_events
set --erase overwrite_hists
set --erase show_plots
