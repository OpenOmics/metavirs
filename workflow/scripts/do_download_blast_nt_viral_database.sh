#!/usr/bin/env bash
set -euo pipefail

function err() { cat <<< "$@" 1>&2; }
function fatal() { cat <<< "$@" 1>&2; exit 1; }
function retry() {
    # Tries to run a cmd 5 times before failing 
    # If a command is successful, it will break out of attempt loop
    # Failed attempts are padding with the following exponential 
    # back-off strategy {4, 16, 64, 256, 1024} in seconds
    # @INPUTS "$@"" = cmd to run 
    # @CALLS fatal() if command cannot be run in 5 attempts
    local n=1
    local max=5
    local attempt=true # flag for while loop
    while $attempt; do
        # Attempt command and break out of attempt loop if successful 
        "$@" && attempt=false || {
        # Try again up to 5 times 
        if [[ $n -le $max ]]; then
            err "Command failed: $@"
            delay=$(( 4**$n ))
            err "Attempt: ${n}/${max}. Trying again in ${delay} seconds!\n"
            sleep $delay;
            ((n++))
        else
            fatal "Fatal: command failed after max retry attempts!"
        fi
      }
    done
}


function extract () {
    # Extract/uncompress different types of files,
    # @INPUT $1 = File to extract.
    if [ -f $1 ] ; then
        case $1 in
            *.tar.bz2)   tar vxjf $1    ;;
            *.tar.gz)    tar vxzf $1    ;;
            *.bz2)       bunzip2 $1     ;;
            *.rar)       rar x $1       ;;
            *.gz)        gunzip $1      ;;
            *.tar)       tar vxf $1     ;;
            *.tbz2)      tar vxjf $1    ;;
            *.tgz)       tar vxzf $1    ;;
            *.zip)       unzip $1       ;;
            *.Z)         uncompress $1  ;;
            *.7z)        7z x $1    ;;
            *)           echo "'$1' cannot be extracted via extract()" ;;
        esac
    else
        err "'$1' is not a valid file"
    fi
}


function download() {
    # This function download a given a URL.
    # @INPUT $1 = URL to the json file
    local download_file="$(basename "${1}")"
    local download_log="$(basename "${1%.*}.log")"
    mkdir -p logs
    retry wget -O "${download_file}" "${1}" -P "${PWD}"  > "logs/${download_log}" 2>&1
}


function main() {
    # Main entry point of script
    # Download JSON containing additional
    # files to download
    download "https://ftp.ncbi.nlm.nih.gov/blast/db/nt_viruses-nucl-metadata.json"

    # Download each viral nt database (0-17),
    # https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html
    # Partially non-redundant nucleotide sequences from 
    # all traditional divisions of GenBank, EMBL, and DDBJ 
    # excluding GSS, STS, PAT, EST, HTG, and WGS.
    while read f; do 
        tarball=$(basename "$f")
        echo "Downloading... $f";
        download "${f}"
        echo "Extracting... ${tarball}"
        extract "$tarball"
    done < \
    <( grep --color=never 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/' \
        nt_viruses-nucl-metadata.json \
            | sed 's/,$/\n/g' \
            | sed "s/^[ \t]*//" \
            | sed '/^[[:space:]]*$/d' \
            | sed 's/"//g'
    )
}


# Call main method
main "$@"