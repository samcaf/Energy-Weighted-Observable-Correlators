#!/usr/local/bin/fish

# Writing
for temp in .1 .2 .3 .4 .5
  for Rjet in 1000.0 1.0
    for rsub in .0 .1
    ./write_tools/write_ewocs -n 10000 -l hadron -p qcd -j 2 -s 1 -R $Rjet -r $rsub --pt_min 50 --pt_max 3000 --frag_temp $temp --verbose 1
    end
  end
end

# Plotting
./plot_tools/scripts/fragtemp_qcd.sh
