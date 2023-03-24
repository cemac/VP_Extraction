#! /bin/bash -l
conda activate DRUID_VP
function usage {
          echo "" 1>&2;
          echo "Usage: $0 -s <<YYYYmmdd>> -e <<YYYYmmdd>> -m <<CVP|QVP>> -f <<file>> [-h] [-v]" 1>&2;
          echo "" 1>&2;
          echo "Required Arguments: " 1>&2;
          echo "  -s : Start date for VP extraction" 1>&2;
          echo "  -e : End date for VP extraction" 1>&2;
          echo "  -m : VP mode, either CVP or QVP (case sensitive)" 1>&2;
          echo "" 1>&2;
          echo "Options: " 1>&2;
          echo "  -f : Path to parameters file (optional, currently only used for CVP extraction)" 1>&2;
          echo "  -h : Show this usage helper" 1>&2;
          echo "  -v : Run programs in verbose mode" 1>&2;
          echo "" 1>&2;
          echo "$2" 1>&2;
          exit $1;
}

verbose="false"

while getopts ":s:e:l:m:f:h" flag; do
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
        f)
            f=${OPTARG}
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

#Get site length for specified mode and params file
site_len=$( ./vp_params.py "${m}" "${f}" )

if [ $? != 0 ]; then
    exit
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
python VP_Main.py $m \$SLURM_ARRAY_TASK_ID $date1_formatted $date2_formatted $site_len $f $verbose
EOF

echo "Created batch script for ${m} job, ${date1_formatted}-${date2_formatted}, ${site_len} sites."
echo "Run this batch with:"
echo "sbatch vp_slurm.sb"
