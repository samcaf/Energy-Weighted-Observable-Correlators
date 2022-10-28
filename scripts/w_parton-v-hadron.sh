#!/usr/local/bin/fish

# ==========================
# Getting arguments
# ==========================
./setup/get_args.sh $argv


# ==========================
# Writing
# ==========================
if $generate_events
for Rjet in 1000.0 1.0;
  for rsub in 0.0 0.1;
    # $alg_prefix'kt' with -s 1, $alg_prefix'ca' with -s 0
    ./write_tools/write_ewocs -n 10000 -l parton -p w -j $alg_prefix'akt' -s $alg_prefix'ca' -R $Rjet -r $rsub --pt_min 50 --pt_max 3000;
    ./write_tools/write_ewocs -n 10000 -l hadron -p w -j $alg_prefix'akt' -s $alg_prefix'ca' -R $Rjet -r $rsub --pt_min 50 --pt_max 3000;
  end;
end;
end


# ==========================
# Plotting
# ==========================
if begin ; $overwrite_hists ; and $show_plots ; end
  ./plot_tools/plot_ewocs --nbins $nbins -n 10000 -l parton -p w -j $alg_prefix'akt' -s $alg_prefix'ca' -R 1000 1 -r 0 .1 --pt_min 50 --pt_max 3000 --plot_type pvh --overwrite_hists true --show_plots true
else if $overwrite_hists
  ./plot_tools/plot_ewocs --nbins $nbins -n 10000 -l parton -p w -j $alg_prefix'akt' -s $alg_prefix'ca' -R 1000 1 -r 0 .1 --pt_min 50 --pt_max 3000 --plot_type pvh --overwrite_hists true --show_plots false
else if $show_plots
  ./plot_tools/plot_ewocs --nbins $nbins -n 10000 -l parton -p w -j $alg_prefix'akt' -s $alg_prefix'ca' -R 1000 1 -r 0 .1 --pt_min 50 --pt_max 3000 --plot_type pvh --overwrite_hists false --show_plots true
end

# ==========================
# Resetting arguments
# ==========================
./setup/clear_args.sh
