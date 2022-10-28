#!/usr/local/bin/fish

# ==========================
# Getting arguments
# ==========================
./setup/get_args.sh $argv


# ==========================
# Writing
# ==========================
if $generate_events
for rsub in 0.0 0.2 0.4 0.6 0.8 1.0 1.2
  ./write_tools/write_ewocs -n 10000 -l parton -p top -j $alg_prefix'antikt' -s $alg_prefix'ca' -R 2.0 -r $rsub --pt_min 50 --pt_max 3000 -E 1000 --s_channel gmZ
end
# Note: need both gamma and Z in s channel
end

# ==========================
# Plotting
# ==========================
if begin ; $overwrite_hists ; and $show_plots ; end
  ./plot_tools/plot_ewocs --nbins $nbins -n 10000 -l parton -p top -j $alg_prefix'antikt' -s $alg_prefix'ca' -R 2.0 -r 0 .2 .4 .6 .8 1 1.2 --pt_min 50 --pt_max 3000 -E 1000 --s_channel gmZ --plot_type sub_rads --overwrite_hists true --show_plots true
else if $overwrite_hists
  ./plot_tools/plot_ewocs --nbins $nbins -n 10000 -l parton -p top -j $alg_prefix'antikt' -s $alg_prefix'ca' -R 2.0 -r 0 .2 .4 .6 .8 1 1.2 --pt_min 50 --pt_max 3000 -E 1000 --s_channel gmZ --plot_type sub_rads --overwrite_hists true --show_plots false
else if $show_plots
  ./plot_tools/plot_ewocs --nbins $nbins -n 10000 -l parton -p top -j $alg_prefix'antikt' -s $alg_prefix'ca' -R 2.0 -r 0 .2 .4 .6 .8 1 1.2 --pt_min 50 --pt_max 3000 -E 1000 --s_channel gmZ --plot_type sub_rads --overwrite_hists false --show_plots true
end

# ==========================
# Resetting arguments
# ==========================
./setup/clear_args.sh
