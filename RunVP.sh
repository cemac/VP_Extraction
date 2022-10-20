#! /bin/bash -l

#Change this parameter if the length of the site list in VP_Main.py changes.
#Current values are: in QVP mode - 7; in CVP mode - 24

QVP_len=$( ./vp_params.py QVP )
CVP_len=$( ./vp_params.py CVP )

function usage {
          echo "" 1>&2;
          echo "Usage: $0 -s <<YYYYmmdd>> -e <<YYYYmmdd>> -m <<CVP|QVP>> [-h] [-v]" 1>&2;
          echo "" 1>&2;
          echo "Required Arguments: " 1>&2;
          echo "  -s : start date for vp extraction" 1>&2;
          echo "  -e : end date for vp extraction" 1>&2;
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
if [ -z "${s}" ] || [ -z "${e}" ] || [ -z "${m}" ]; then
    usage 1 "-s -e and -m arguments are all required"
fi

if [ ${m} == "CVP" ]; then
    site_len=$CVP_len
elif [ ${m} == "QVP" ]; then
    site_len=$QVP_len
fi

date1=$( date -d $s +%s )
date2=$( date -d $e +%s )

if  [ $date2 -lt $date1 ]; then
  echo "start date occurs after end date. Swapping dates"
  date2=$( date -d $s +%s )
  date1=$( date -d $e +%s )
fi

date_len=$(( ($date2 - $date1 )/(60*60*24)+1))
date1_formatted=$( date -u -d @${date1} +'%Y%m%d')
date2_formatted=$( date -u -d @${date2} +'%Y%m%d')

Max_iter=$(( $date_len*$site_len ))
mkdir -p Output
cat > vp_slurm.sb <<-EOF
#!/bin/bash -l
#SBATCH --job-name=${m}_job_array    # Job name
#SBATCH --time=12:00:00             # Time limit per array task hrs:min:sec
#SBATCH --output=Output/%j-%A_%a.out       # Standard output and error log

#SBATCH --array=1-$Max_iter              # Array range
source activate DRUID_VP
python VP_Main.py $m \$SLURM_ARRAY_TASK_ID $date1_formatted $date2_formatted $site_len $verbose
EOF

echo "sbatch vp_slurm.sb"
