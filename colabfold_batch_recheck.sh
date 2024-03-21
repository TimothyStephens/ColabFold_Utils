#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
start=`date +%s.%N`


#### Function to help keep track of progress.
function run_cmd(){
        local CMD="$@"
        echo "[`date`]      CMD: $CMD"
        eval $CMD
}
function log(){
        echo "[`date`]      LOG: $@"
}


#### Envs
OUTDIR=""


#### Useage information
usage() {
echo -e "##
## $(basename ${0})
##

Wrapper script to run ColabFold Batch on Search results and check that the results were correctly generated.

Sometimes colabfold_batch will fail with exitstatus==0 if at startup TensorFlow fails to load,
which makes it challenging to know if all your jobs finished correctly. This wrapper checks for
the presence of the correct output files, which should make this analysis more robust and scriptable.

Usage: 
$(basename $0) --OUTDIR_TO_CHECK results [colabfold_batch options] input results

Options:
--OUTDIR_TO_CHECK          The 'results' directory passed to colabfold_batch to check for completness (Required)

-h, --help                 This help message
--debug                    Run debug mode

##
## colabfold_batch help message:
##
" 1>&2
}


# See https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
set +eu
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    --OUTDIR_TO_CHECK)
      OUTDIR="$2"
      shift # past argument
      shift # past value
      ;;
    -h|--help)
      usage
      POSITIONAL+=("$1") # save it in an array for later
      shift # past argument
      ;;
    --debug)
      set -x
      shift # past argument
      ;;
    *) # unknown option
      POSITIONAL+=("$1") # save it in an array for later
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters
set -eu


#### Check that $OUTDIR is set
if [[ "$OUTDIR" == "" ]];
then
  log "COMPLETNESS CHECK [ERROR]: '--OUTDIR_TO_CHECK' was not set! Please make this value the same as colabfold_batch results dir."
  exit 1
else
  log "COMPLETNESS CHECK [OK]: Will be checking for the correct output of colabfold_batch in ${OUTDIR}"
fi

#### Run colabfold_search
run_cmd "colabfold_batch $@"


#### Check if colabfold_batch finished correctly
# Check if output dir was created (wont exist if colabfold_batch failed during library loading)
if [[ ! -d "${OUTDIR}" ]];
then
  log "COMPLETNESS CHECK [ERROR]: ${OUTDIR} does not exist! colabfold_batch might not have finished correctly."
  exit 1
else
  log "COMPLETNESS CHECK [OK]: ${OUTDIR} exist!"
fi

## Get *.a3m files in output results dir and check if each has a *.done.txt file
##  - The presence of the *.done.txt file signifies that script reached the end of its run.
# Check if any *.a3m files were found
if [[ $(find "${OUTDIR}" -name "*.a3m" | wc -l) -eq 0 ]];
then
  log "COMPLETNESS CHECK [ERROR]: No *.a3m files found in ${OUTDIR}! colabfold_batch might not have finished correctly."
  exit 1
else
  log "COMPLETNESS CHECK [OK]: *.a3m files found in ${OUTDIR}!"
fi
# For each *.a3m file, check if *.done.txt file exists
for A3M in $(find "${OUTDIR}" -name "*.a3m" | sort);
do
  if [[ ! -e "${A3M%*.a3m}.done.txt" ]];
  then
    log "COMPLETNESS CHECK [ERROR]: No *.done.txt file found for ${A3M}! colabfold_batch might not have finished correctly."
    exit 1
  fi
done
log "COMPLETNESS CHECK [OK]: All *.a3m files found in ${OUTDIR} have *.done.txt files!"

## Looks like everything passed.
log "Looks like colabfold_batch created all the expected output files"


#### Print run time
end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
dhms=$(printf '%02dd:%02dh:%02dm:%02fs\n' $(echo -e "$runtime/86400\n$runtime%86400/3600\n$runtime%3600/60\n$runtime%60"| bc))
log "Finished running ColabFold Batch! It ran for: ${dhms}"
