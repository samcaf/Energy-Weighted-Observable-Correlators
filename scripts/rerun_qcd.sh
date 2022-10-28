#!/usr/local/bin/fish

# =============================
# QCD
# =============================
echo -en "Start time : "(date +%m/%d-%H:%M)"\n\n"

# -------------------------
# Subjet radii
# -------------------------
# Parton
./scripts/qcd_parton_subr.sh $argv
./scripts/qcd_parton_subr_10k.sh $argv
# Ecm = 10TeV (rather than the default of 4TeV)

# Hadron
./scripts/qcd_hadron_subr.sh $argv
# No s-channel Z boson
./scripts/qcd_parton_noZ_subr_10k.sh $argv
./scripts/qcd_parton_noZ_subr_10k_kt.sh $argv

# -------------------------
# Parton vs. Hadron
# -------------------------
./scripts/qcd_parton-v-hadron.sh $argv
./scripts/qcd_parton-v-hadron_10k.sh $argv

# -------------------------
# Fragmentation Temperature
# -------------------------
./scripts/qcd_fragtemp.sh $argv
./scripts/qcd_fragtemp_r0.sh $argv

echo -en "\n\nEnd time : "(date +%m/%d-%H:%M)
