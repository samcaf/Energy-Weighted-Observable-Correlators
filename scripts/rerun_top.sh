#!/usr/local/bin/fish

# =============================
# Top
# =============================
echo -en "Start time : "(date +%m/%d-%H:%M)"\n\n"

# -------------------------
# Subjet radii
# -------------------------
./scripts/top_parton_subr.sh $argv
./scripts/top_parton_subr_500.sh $argv
./scripts/top_parton_subr_1000.sh $argv

echo -en "\n\nEnd time : "(date +%m/%d-%H:%M)
