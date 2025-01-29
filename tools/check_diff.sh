#!/usr/bin/env bash

# Parameters
inputFile=$1
outputFile=$2
mode=$3
commandPath=$4

# Initial check
if test -f "${outputFile}"; then
  # Destination file exists
  if cmp -s ${inputFile} ${outputFile}; then
    # Destination file and source file are similar, exiting
    exit 0
  fi
fi

# Prepare destination directory
dstDir=$(dirname "${outputFile}")
mkdir -p ${dstDir}

# Copy file
cp -f ${inputFile} ${outputFile}.tmp

# Get file extension
srcFile=$(basename "${inputFile}")
srcExt="${srcFile##*.}"

if test "${srcExt}" = "h" -o "${srcExt}" = "cc"; then
  if test "${mode}" = "ECSABER"; then
    # ECSABER update
    if test "${srcFile}" = "Variables.h"; then
      sed -i -e s/"Variables"/"JediVariables"/g ${outputFile}.tmp
      sed -i -e s/"JediVariablesBase"/"VariablesBase"/g ${outputFile}.tmp
    elif test "${srcFile}" = "Variables.cc"; then
      sed -i -e s/"Variables"/"JediVariables"/g ${outputFile}.tmp
      sed -i -e s/"JediVariablesBase"/"VariablesBase"/g ${outputFile}.tmp
      sed -i -e s/"\/JediVariables"/"\/Variables"/g ${outputFile}.tmp
    else
      sed -i -e s/"Variables::Variables"/"JediVariables::JediVariables"/g ${outputFile}.tmp
      sed -i -e s/"oops::Variables"/"oops::JediVariables"/g ${outputFile}.tmp
      sed -i -e s/" Variables"/" JediVariables"/g ${outputFile}.tmp
      sed -i -e s/"<Variables>"/"<JediVariables>"/g ${outputFile}.tmp
      sed -i -e s/"^Variables"/"JediVariables"/g ${outputFile}.tmp
    fi
  elif test "${mode}" = "QUENCHXX"; then
    # QUENCH update
    sed -i -e s/"quench"/"quenchxx"/g ${outputFile}.tmp
    sed -i -e s/"src\/"/"quenchxx\/"/g ${outputFile}.tmp
  else
    echo "Wrong mode in check_diff.sh: "${mode}
    exit 1
  fi
fi

# Check patch size
if test -f "${outputFile}.patch"; then
  patchLength=`cat ${outputFile}.patch | wc -l`
  if test "${patchLength}" = "0"; then
    rm -f ${outputFile}.patch
  fi
fi

if test -f "${outputFile}.patch"; then
  # Apply residual patch
  cp -f ${outputFile}.tmp ${outputFile}.tmp.bak
  patch -s ${outputFile}.tmp ${outputFile}.patch

  # Compare and update if needed
  if cmp -s ${outputFile}.tmp ${outputFile}; then
    rm -f ${outputFile}.tmp ${outputFile}.tmp.bak
  else
    echo "--  - Update needed for: ${outputFile}"
    echo "echo -e \"--  - \e[31m[with patch] Difference for ${outputFile}: [u]pdate, [m]eld or keep (any key)?\e[0m\"; read action; if test \"\${action}\" = \"u\"; then cp -f ${outputFile}.tmp ${outputFile}; elif test \"\${action}\" = \"m\"; then meld ${outputFile}.tmp ${outputFile} ${outputFile}.tmp.bak; fi; diff -u ${outputFile}.tmp.bak ${outputFile} > ${outputFile}.patch; patchLength=\`cat ${outputFile}.patch | wc -l\`; if test "\$\{patchLength\}" = "0"; then rm -f ${outputFile}.patch;fi; rm -f ${outputFile}.tmp ${outputFile}.tmp.bak" >> ${commandPath}
    echo "" >> ${commandPath}
 fi
else
  if test -f "${outputFile}"; then
    # Create patch if needed
    diff -u ${outputFile}.tmp ${outputFile} > ${outputFile}.patch
    patchLength=`cat ${outputFile}.patch | wc -l`
    if test "${patchLength}" = "0"; then 
      rm -f ${outputFile}.patch
      rm -f ${outputFile}.tmp
    else
      echo "--  - New patch needed for: "${outputFile}
      rm -f ${outputFile}.patch
      echo "echo -e \"--  - \e[92m[no patch] Difference for ${outputFile}: [u]pdate, [m]eld or keep (any key)?\e[0m\"; read action; if test \"\${action}\" = \"u\"; then cp -f ${outputFile}.tmp ${outputFile}; elif test \"\${action}\" = \"m\"; then meld ${outputFile}.tmp ${outputFile}; fi; diff -u ${outputFile}.tmp ${outputFile} > ${outputFile}.patch; patchLength=\`cat ${outputFile}.patch | wc -l\`; if test "\$\{patchLength\}" = "0"; then rm -f ${outputFile}.patch;fi;rm -f ${outputFile}.tmp" >> ${commandPath}
      echo "" >> ${commandPath}
    fi
  else
    # New file
    mv ${outputFile}.tmp ${outputFile}
  fi
fi
