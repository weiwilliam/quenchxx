#!/usr/bin/env bash

# Parameters
input=$1
output=$2
commandPath=$3

# Script directory
scriptDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Loop recursively over items
for item in `find ${input}`; do
  # Check if item is a file
  if test -f "${item}"; then
    # Get input directory and file
    dir=$(dirname "${item}")
    file=$(basename "${item}")

    # Update file
    bash ${scriptDir}/check_diff.sh ${item} ${output}/${file} ${commandPath}
  fi
done
