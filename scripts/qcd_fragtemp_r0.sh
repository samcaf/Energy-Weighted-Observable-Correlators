#!/usr/local/bin/fish

# ==========================
# Getting arguments
# ==========================
./setup/get_args.sh $argv


# ==========================
# Writing
# ==========================
if $generate_events
for temp in .1 .2 .3 .4 .5
  ./write_tools/write_ewocs -n 10000 -l hadron -p qcd -j $alg_prefix'akt' -s $alg_prefix'ca' -R 1.0 -r 0.0 --pt_min 50 --pt_max 3000 --frag_temp $temp --verbose 2
end
end


# ==========================
# Plotting
# ==========================
if begin ; $overwrite_hists ; and $show_plots ; end
./plot_tools/plot_ewocs --nbins $nbins -n 10000 -l hadron -p qcd -j $alg_prefix'akt' -s $alg_prefix'ca' -R 1.0 -r 0.0 --pt_min 50 --pt_max 3000 --frag_temp .1 .2 .3 .4 .5 --plot_type thermal --overwrite_hists true --show_plots true
else if $overwrite_hists
./plot_tools/plot_ewocs --nbins $nbins -n 10000 -l hadron -p qcd -j $alg_prefix'akt' -s $alg_prefix'ca' -R 1.0 -r 0.0 --pt_min 50 --pt_max 3000 --frag_temp .1 .2 .3 .4 .5 --plot_type thermal --overwrite_hists true --show_plots false
else if $show_plots
./plot_tools/plot_ewocs --nbins $nbins -n 10000 -l hadron -p qcd -j $alg_prefix'akt' -s $alg_prefix'ca' -R 1.0 -r 0.0 --pt_min 50 --pt_max 3000 --frag_temp .1 .2 .3 .4 .5 --plot_type thermal --overwrite_hists false --show_plots true
end


# ==========================
# Resetting arguments
# ==========================
./setup/clear_args.sh
