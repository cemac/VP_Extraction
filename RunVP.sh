#! /bin/bash -l

#Change this parameter if the length of the site list in VP_Main.py changes.
#Current values are: in QVP mode - 7; in CVP mode - 24

QVP_len=$( ./vp_params.py QVP )
CVP_len=$( ./vp_params.py CVP )

if [ $QVP_len -gt $CVP_len ]; then
  max_len=$QVP_len
else
  max_len=$CVP_len
fi

function usage {
          echo "" 1>&2;
          echo "Usage: $0 -s <<YYYYmmdd>> -e <<YYYYmmdd>> -l << int >> -m <<CVP|QVP>> [-h] [-v]" 1>&2;
          echo "" 1>&2;
          echo "Required Arguments: " 1>&2;
          echo "  -s : start date for vp extraction" 1>&2;
          echo "  -e : end date for vp extraction" 1>&2;
          echo "  -l : length of site list, minimum 1 maximum $QVP_len in QVP mode, or $CVP_len in CVP mode. List can be editted in VP_Main.py" 1>&2;
          echo "  -m : VP mode, either CVP or QVP (case sensitive)" 1>&2;
          echo "" 1>&2;
          echo "Options: " 1>&2;
          echo "  -h : Show this usage helper" 1>&2;
          echo "  -v : Run programs in verbose mode" 1>&2;
          echo "" 1>&2;
          echo "$2" 1>&2;
          exit $1;
}

verbose="false"

while getopts ":s:e:l:m:h" flag; do
    case "${flag}" in
        l)
            l=${OPTARG}
            ((l >= 1 && l <=$max_len )) || usage 1 "site list length outside of accepted bounds"
            ;;
        s)
            s=${OPTARG}
            ;;
        e)
            e=${OPTARG}
            ;;
        m)
            m=${OPTARG}
            ((m == 'CVP' || m == 'QVP')) || usage 1 "VP mode must be either QVP or CVP"
            ;;
        h)
            usage 0 ""
            ;;
        v)
            verbose="true"
            ;;
        *)
            usage 1 "Unrecognised option ${flag}"
            ;;
    esac
done
shift $((OPTIND-1))

#Requires all three arguments
if [ -z "${s}" ] || [ -z "${e}" ] || [ -z "${l}" ] || [ -z "${m}" ]; then
    usage 1 "-s -e -l and -m arguments are all required"
fi

if [ ${m} == "CVP" ] &&  [ ${l} -gt $CVP_len ]; then
    usage 1 "Maximum value of site length in CVP mode is $CVP_len"
elif [ ${m} == "QVP" ] && [ ${l} -gt $QVP_len ]; then
    usage 1 "Maximum value of site length in QVP mode is $QVP_len"
fi

date1=$( date -d $s +%s )
date2=$( date -d $e +%s )

if  [ $date2 -lt $date1 ]; then
  echo "start date occurs after end date. Swapping dates"
  date2=$( date -d $s +%s )
  date1=$( date -d $e +%s )
fi

date_len=$(( ($date2 - $date1 )/(60*60*24) + 1 ))
date1_formatted=$( date -u -d @${date1} +"%Y%m%d" )
date2_formatted=$( date -u -d @${date2} +"%Y%m%d" )

Max_iter=$(( $date_len*$l ))

cat > vp_slurm.sb <<-EOF
#!/bin/bash -l
#SBATCH --job-name=${m}_job_array    # Job name
#SBATCH --time=12:00:00             # Time limit per array task hrs:min:sec
#SBATCH --output=%j-%A_%a.out       # Standard output and error log
#SBATCH --array=1-$Max_iter              # Array range

source activate DRUID_VP
python VP_Main.py $m \$SLURM_ARRAY_TASK_ID $date1_formatted $date2_formatted $l $verbose
EOF

echo "sbatch vp_slurm.sb"
