#!/usr/local/bin/fish

set qcd_cmd "./scripts/rerun_qcd.sh $argv"
set w_cmd "./scripts/rerun_w.sh $argv"
set top_cmd "./scripts/rerun_top.sh $argv"

# MacOS Scripts
osascript -e "tell application \"Terminal\" to do script \"cd $PWD; $qcd_cmd\""
osascript -e "tell application \"Terminal\" to do script \"cd $PWD; $w_cmd\""
osascript -e "tell application \"Terminal\" to do script \"cd $PWD; $top_cmd\""
