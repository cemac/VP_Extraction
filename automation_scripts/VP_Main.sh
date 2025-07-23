#! /bin/bash -l
r=$1
s=$2
c=$3
index=$4
verbose=$5

args=''
if [ verbose ]; then args='-v'; fi
echo 'VP_Main.sh' $r $s $c $index $args
this_date=$s
if (( $index > 0 )); then this_date=$(date -d "$s + $index day" +'%Y%m%d'); fi
echo 'about to run vp_extraction' $r $this_date
python vp_extraction.py -r $r -t $this_date -c $c $args
