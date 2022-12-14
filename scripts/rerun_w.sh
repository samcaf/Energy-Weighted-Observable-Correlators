#!/usr/local/bin/fish


# =============================
# W
# =============================
echo -en "Start time : "(date +%m/%d-%H:%M)"\n\n"

# -------------------------
# Subjet radii
# -------------------------
./scripts/w_parton_subr.sh $argv
./scripts/w_parton_subr_500.sh $argv
# Again, suffix is Ecm in GeV

# -------------------------
# Parton vs. Hadron
# -------------------------
./scripts/w_parton-v-hadron.sh $argv
./scripts/w_parton-v-hadron_500.sh $argv

echo -en "\n\nEnd time : "(date +%m/%d-%H:%M)
