#!/usr/local/bin/fish

# Writing
for Rjet in 1000.0 1.0;
  for rsub in 0.0 0.1;
    # kt with -s 1, ca with -s 0
    ./write_tools/write_ewocs -n 10000 -l parton -p w -j akt -s ca -R $Rjet -r $rsub --pt_min 50 --pt_max 3000;
    ./write_tools/write_ewocs -n 10000 -l hadron -p w -j akt -s ca -R $Rjet -r $rsub --pt_min 50 --pt_max 3000;
  end;
end;

# Plotting
./plot_tools/plot_ewocs -n 10000 -l parton -p w -j akt -s ca -R 1000 1 -r 0 .1 --pt_min 50 --pt_max 3000 --plot_type pvh
