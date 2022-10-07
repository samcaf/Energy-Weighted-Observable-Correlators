#!/usr/local/bin/fish

for temp in .1 .2 .3 .4 .5
  ./pythia_ewocs -n 10000 -l hadron -p qcd -j 2 -s 1 -R 1.0 -r 0.0 --pt_min 50 --pt_max 3000 --frag_temp $temp --verbose 2
end

cd python; ./scripts/fragtemp_qcd_r0.sh
cd ..
