#!/usr/local/bin/fish

# Writing
for temp in .1 .2 .3 .4 .5
  ./write_tools/write_ewocs -n 10000 -l hadron -p qcd -j akt -s ca -R 1.0 -r 0.0 --pt_min 50 --pt_max 3000 --frag_temp $temp --verbose 2
end

# Plotting
./plot_tools/plot_ewocs -n 10000 -l hadron -p qcd -j akt -s ca -R 1.0 -r 0.0 --pt_min 50 --pt_max 3000 --frag_temp .1 .2 .3 .4 .5 --plot_type thermal
