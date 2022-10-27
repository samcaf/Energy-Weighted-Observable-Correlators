#!/usr/local/bin/fish

# Resetting scripts and erasing global variables:
set --erase ee_algs
find scripts/ -type f -print0 | xargs -0 sed -i '' -E 's|' ee_kt'|' kt'|g'
find scripts/ -type f -print0 | xargs -0 sed -i '' -E 's|' ee_akt'|' akt'|g'
find scripts/ -type f -print0 | xargs -0 sed -i '' -E 's|' ee_ca'|' ca'|g'
set --erase generate_events
set --erase overwrite_hists
set --erase show_plots
