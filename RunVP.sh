#! /bin/bash -l
conda activate DRUID_VP
function usage {
          echo "" 1>&2;
          echo "Usage: $0 -s <<YYYYmmdd>> -e <<YYYYmmdd>> -m <<CVP|QVP>> -c <<cfg_file>> [-h] [-v]" 1>&2;
          echo "" 1>&2;
          echo "Required Arguments: " 1>&2;
          echo "  -r : radar name" 1>&2;
          echo "  -s : Start date for VP extraction" 1>&2;
          echo "  -e : End date for VP extraction" 1>&2;
          echo "  -c : Path to config file" 1>&2;
          echo "" 1>&2;
          echo "Options: " 1>&2;
          echo "  -h : Show this usage helper" 1>&2;
          echo "  -v : Run programs in verbose mode" 1>&2;
          echo "" 1>&2;
          echo "$2" 1>&2;
          exit $1;
}

verbose="false"

while getopts ":r:s:e:c" flag; do
    case "${flag}" in
        r)
            r=${OPTARG}
            echo 'r='$r
            ;;
        s)
            s=${OPTARG}
            echo 's='$s
            ;;
        e)
            e=${OPTARG}
            echo 'e='$e
            ;;
        c)
            c=${OPTARG}
            echo 'c='$c
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

#Requires all these arguments
if [ -z "${r}" ] || [ -z "${s}" ] || [ -z "${e}" ] ||  [ -z "${c}" ]; then
    usage 1 "-r -s -e and -c arguments are all required"
fi

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
#date2_formatted=$( date -u -d @${date2} +'%Y%m%d')

Max_iter=$(( $date_len ))
mkdir -p Output
this_date=$date1_formatted
for ((i=1; i<=Max_Iter; i++))
do
    args=''
    if [ verbose ]; then args='-v'; fi
    
    sbatch --account=ncas_radar --partition=standard --time=04:00:00 --output=Output/$r_${this_date}.out  --job-name=$r_${this_date} --wrap="vp_extraction.py -r $r -t $this_date -c $c $args" 
    this_date=$(date +"%Y%m%d" -d "$this_date + 1 day")

done


