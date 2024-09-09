#!/bin/bash
base_dir=$PWD
 
find "$base_dir" -name 'ABCstat.txt' -exec sh -c '
  for file do
    line_count=$(wc -l < "$file")
    if [ "$line_count" -ne 500 ]; then
      echo $(dirname "$file") | grep topo
    fi
  done
' sh {} + > unfinished_sim.txt

