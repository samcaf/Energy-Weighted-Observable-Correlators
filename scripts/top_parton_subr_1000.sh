#!/usr/local/bin/fish

# Writing
for rsub in 0.0 0.2 0.4 0.6 0.8 1.0 1.2
  ./write_tools/write_ewocs -n 10000 -l parton -p top -j antikt -s ca -R 2.0 -r $rsub --pt_min 50 --pt_max 3000 -E 1000 --s_channel gmZ
end
# Note: need both gamma and Z in s channel

# Plotting
./plot_tools/scripts/subrads_top_parton_1000.sh
