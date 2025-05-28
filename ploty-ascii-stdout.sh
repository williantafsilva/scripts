#!/bin/bash
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Print a simple Y plot using ASCII.
##Requirements: gnuplot.

##Input $1: Y data (comma-separated list of numbers; use | paste -sd,)
##Output: ASCII plot.

##Usage: 
##ploty-ascii-stdout.sh <Y DATA>

############################################################################
##ACTIONS:

##Input.

RUNDATE=$(date +"%Y%m%d%H%M%S")
YVALUES=($(echo $1 | sed 's/,/\n/g')) ##Convert comma-separated list to array.

##Process.

##Load modules.
module load gnuplot

##Create a temporary data file with indices and values (x, y).
TMPFILE=$(echo "tmp-ploty-ascii-${RUNDATE}.txt")
INDEX=1
for VALUE in "${YVALUES[@]}" ; do
  echo "${INDEX} ${VALUE}" >> "${TMPFILE}"
  ((INDEX++))
done

##Plot the data using Gnuplot.
gnuplot -persist << EOF
set terminal dumb
set title "Index X Value"
set xlabel "Index"
set ylabel "Value"
#set grid

#Determine axis range with some padding.
stats "${TMPFILE}" using 1 nooutput
x_min = STATS_min - 0.5
x_max = STATS_max + 0.5

stats "${TMPFILE}" using 2 nooutput
y_min = STATS_min - (STATS_max - STATS_min) * 0.1
y_max = STATS_max + (STATS_max - STATS_min) * 0.1

set xrange [x_min:x_max]
set yrange [y_min:y_max]

plot "${TMPFILE}" using 1:2 with points pointtype 15 notitle
EOF

##Remove the temporary file.
rm "${TMPFILE}"
