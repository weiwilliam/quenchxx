#!/usr/bin/env bash

# Parameters
checkDiff=$1
inputDir=$2
outputDir=$3
mode=$4
commandPath=$5

# Loop over items
for item in `find ${inputDir}`; do
  # Check if item is a file
  if test -f "${item}"; then
    # Get full directories and file name
    fullInputDir=$(dirname "${item}")
    fullOutputDir=${fullInputDir/"${inputDir}"/"${outputDir}"}
    file=$(basename "${item}")

    # Update file
    bash ${checkDiff} ${item} ${fullOutputDir}/${file} ${mode} ${commandPath}
  fi
done
