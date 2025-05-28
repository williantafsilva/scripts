#!/bin/bash
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Print a simple X,Y plot using ASCII.
##Requirements: gnuplot.

##Input $1: X label.
##Input $2: X data (comma-separated list of numbers; use | paste -sd,)
##Input $3: Y label.
##Input $4: Y data (comma-separated list of numbers; use | paste -sd,)
##Output: ASCII plot.

##Usage: 
##plotxy-ascii-stdout.sh <X LABEL> <X DATA> <Y LABEL> <Y DATA>

############################################################################
##ACTIONS:

##Input.

RUNDATE=$(date +"%Y%m%d%H%M%S")
XLABEL="$1"
XVALUES=($(echo $2 | sed 's/,/\n/g')) ##Convert comma-separated list to array.
YLABEL="$3"
YVALUES=($(echo $4 | sed 's/,/\n/g')) ##Convert comma-separated list to array.

##Process.

##Load modules.
module load gnuplot

##Create a temporary data file.
TMPFILE=$(echo "tmp-plotxy-ascii-${RUNDATE}.txt")
for i in "${!XVALUES[@]}" ; do
    echo "${XVALUES[$i]} ${YVALUES[$i]}" >> "${TMPFILE}"
done

##Plot the data using Gnuplot.
gnuplot -persist << EOF
set terminal dumb
set title "${XLABEL} X ${YLABEL}"
set xlabel "${XLABEL}"
set ylabel "${YLABEL}"
#set grid

#Determine axis range with some padding.
stats "${TMPFILE}" using 1 nooutput
x_min = STATS_min - (STATS_max - STATS_min) * 0.1
x_max = STATS_max + (STATS_max - STATS_min) * 0.1

stats "${TMPFILE}" using 2 nooutput
y_min = STATS_min - (STATS_max - STATS_min) * 0.1
y_max = STATS_max + (STATS_max - STATS_min) * 0.1

set xrange [x_min:x_max]
set yrange [y_min:y_max]

plot "${TMPFILE}" using 1:2 with points pointtype 15 notitle
EOF

##Remove the temporary file.
rm "${TMPFILE}"
