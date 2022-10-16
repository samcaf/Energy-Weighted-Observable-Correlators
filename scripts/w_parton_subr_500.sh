#!/usr/local/bin/fish

# Writing
for rsub in 0.0 0.1 0.2 0.3 0.4 0.5
  ./write_tools/write_ewocs -n 10000 -l parton -p w -j akt -s ca -R 1.0 -r $rsub --pt_min 50 --pt_max 3000 -E 500
end

# Plotting
./plot_tools/plot_ewocs -n 10000 -l parton -p w -j akt -s ca -R 1.0 -r 0 .1 .2 .3 .4 .5 --pt_min 50 --pt_max 3000 -E 500 --plot_type sub_rads
