#!/usr/local/bin/fish

for Rjet in 1000.0 1.0;
  for rsub in 0.0 0.1;
    # kt with -s 1, ca with -s 0
    ./pythia_ewocs -n 10000 -l parton -p qcd -j 2 -s 1 -R $Rjet -r $rsub --pt_min 50 --pt_max 3000 -E 10000;
    ./pythia_ewocs -n 10000 -l hadron -p qcd -j 2 -s 1 -R $Rjet -r $rsub --pt_min 50 --pt_max 3000 -E 10000;
  end;
end;

cd python; ./scripts/parton-v-hadron_qcd_10k.sh;
cd ..
