#!/usr/local/bin/fish

for rsub in 0.0 0.1 0.2 0.3 0.4 0.5
  ./pythia_ewocs -n 10000 -l parton -p qcd -j 2 -s 1 -R 1.0 -r $rsub --pt_min 50 --pt_max 3000 -E 10000
end

cd python; ./scripts/subrads_qcd_parton_10k.sh
cd ..
