#!/usr/bin/env bash

## Overview: this script produces the index in this directory.


ThisDir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# The name of the index file. If it exists already, delete it.
Index="$ThisDir"/'README.txt'
if [ -f "$Index" ]; then
  rm "$Index"
fi

for file in "$ThisDir"/*; do

  # Skip any file whose name ends in '~' (automatic backup files in Ubuntu) or
  # '.pyc' (compiled python files).
  if [ "${file: -1}" == "~" ] || [ "${file: -4}" == '.pyc' ]; then
    continue
  fi

  # Skip the index file itself
  if [ "$file" == "$Index" ]; then
    continue
  fi

  # Check that a line begins '## Overview:'
  if ! grep -e '^## Overview:' "$file" > /dev/null; then
    echo `basename "$file"` 'is missing a line beginning "## Overview:"' 1>&2
    continue
  fi

  # Print the file name into the index
  echo -e `basename "$file"` "\n" >> "$Index"

  # Print from the line starting '## Overview:' until either a blank line or a
  # line containing only two hashes.
  awk '/^## Overview:/ , /^$/ || /^##$/' "$file" | sed '$d' >> "$Index"
  echo -e "\n\n" >> "$Index"

done
