#!/usr/bin/env bash

# Source and destination paths
srcPath=$1
dstPath=$2
commandPath=$3

# Initial check
if test -f "${dstPath}"; then
  # Destination file exists
  if cmp -s ${srcPath} ${dstPath}; then
    # Destination file and source file are similar, exiting
    exit 0
  fi
fi

# Prepare destination directory
dstDir=$(dirname "${dstPath}")
mkdir -p ${dstDir}

# Copy file
cp -f ${srcPath} ${dstPath}.tmp

# Get file extension
srcFile=$(basename "${srcPath}")
srcExt="${srcFile##*.}"


if test "${srcExt}" = "h" -o "${srcExt}" = "cc"; then
  # Update content
  sed -i -e s/"quench"/"quenchxx"/g ${dstPath}.tmp
  sed -i -e s/"src\/"/"quenchxx\/"/g ${dstPath}.tmp
fi

# Check patch size
if test -f "${dstPath}.patch"; then
  patchLength=`cat ${dstPath}.patch | wc -l`
  if test "${patchLength}" = "0"; then
    rm -f ${dstPath}.patch
  fi
fi

if test -f "${dstPath}.patch"; then
  # Apply residual patch
  cp -f ${dstPath}.tmp ${dstPath}.tmp.bak
  patch -s ${dstPath}.tmp ${dstPath}.patch

  # Compare and update if needed
  if cmp -s ${dstPath}.tmp ${dstPath}; then
    rm -f ${dstPath}.tmp ${dstPath}.tmp.bak
  else
    echo "--  - Update needed for: ${dstPath}"
    echo "echo -e \"--  - \e[31m[with patch] Difference for ${dstPath}: [u]pdate, [m]eld or keep (any key)?\e[0m\"; read action; if test \"\${action}\" = \"u\"; then cp -f ${dstPath}.tmp ${dstPath}; elif test \"\${action}\" = \"m\"; then meld ${dstPath}.tmp ${dstPath} ${dstPath}.tmp.bak; fi; diff -u ${dstPath}.tmp.bak ${dstPath} > ${dstPath}.patch; patchLength=\`cat ${dstPath}.patch | wc -l\`; if test "\$\{patchLength\}" = "0"; then rm -f ${dstPath}.patch;fi; rm -f ${dstPath}.tmp ${dstPath}.tmp.bak" >> ${commandPath}
    echo "" >> ${commandPath}
 fi
else
  if test -f "${dstPath}"; then
    # Create patch if needed
    diff -u ${dstPath}.tmp ${dstPath} > ${dstPath}.patch
    patchLength=`cat ${dstPath}.patch | wc -l`
    if test "${patchLength}" = "0"; then 
      rm -f ${dstPath}.patch
      rm -f ${dstPath}.tmp
    else
      echo "--  - New patch needed for: "${dstPath}
      rm -f ${dstPath}.patch
      echo "echo -e \"--  - \e[92m[no patch] Difference for ${dstPath}: [u]pdate, [m]eld or keep (any key)?\e[0m\"; read action; if test \"\${action}\" = \"u\"; then cp -f ${dstPath}.tmp ${dstPath}; elif test \"\${action}\" = \"m\"; then meld ${dstPath}.tmp ${dstPath}; fi; diff -u ${dstPath}.tmp ${dstPath} > ${dstPath}.patch; patchLength=\`cat ${dstPath}.patch | wc -l\`; if test "\$\{patchLength\}" = "0"; then rm -f ${dstPath}.patch;fi;rm -f ${dstPath}.tmp" >> ${commandPath}
      echo "" >> ${commandPath}
    fi
  else
    # New file
    mv ${dstPath}.tmp ${dstPath}
  fi
fi
