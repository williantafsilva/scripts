#!/bin/bash
############################################################################
############################ STDOUT FUNCTION ###############################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Show useful commands.

##Input: NONE
##Output: List of useful commands.

##Usage: 
##help-stdout.sh 

############################################################################
##ACTIONS:

echo '----------------------------> USEFUL COMMANDS <-----------------------------
---> Full path to file/directory: 
readlink -f <FILE/DIRECTORY>

---> Count lines:
<INPUT> | wc -l

---> List all variables and functions:
declare | less

---> BASH script arguments:
Script location: $0
First argument: $1
Second argument: $2
Last argument: "${@:$#}" or "${@: -1}"
All but last argument: "${@:1:$#-1}"
From third argument to last: "${@:3}"
All arguments: "$@"

---> Screen sessions:
List available screens: screen -ls
Create a screen session: screen -D -R -S <SCREEN NAME>
Attach screen session: screen -r <SESSION ID>
Dettach screen session: screen -d -S <SESSION ID>
Quit screen session: screen -XS <SESSION ID> quit
Quit current screen session: <CTRL+SHIFT+a+d>

---> Extract information from path to file ():
File location: ${<PATHTOFILE>%/*}
File name: ${<PATHTOFILE>##*/}

---> Useful sed commands:
Delete everything before and with pattern: <INPUT> | sed -e '"'"'s/.*PATTERN//'"'"'
Delete everything before pattern but not the pattern: <INPUT> | sed '"'"'s/.*\(PATTERN.*$\)/\1/'"'"'

---> Select lines than match a pattern in columns 4 and 5, and print values in column 3 of the selected lines.
<INPUT> | awk -v L="${L}" '"'"'($4 == L && $5 == "R2") {print $3}'"'"'

---> Select line by number (15th line).
<INPUT> | awk '"'"'NR==15'"'"'

---> Select lines by condition in column 4.
<INPUT> | awk '"'"'$4 < 0.01'"'"'

---> Sort by columns (2nd column alphanumeric, 3rd column numeric, 4th column numeric).
<INPUT> | sort -V -k2,2 -k3,3n -k4,4n

---> Select lines based on condition (tab-separated, define awk variable S, 1st column equals variable S and 5th column equals "R1").
<INPUT> | awk -F '"'"'\t'"'"' -v S="${MYVARIABLE}" '"'"'$1 == S && $5 == "R1"'"'"'

---> Select line based on the presence of a variable followed by tab.
<INPUT> | grep -P "${MYVARIABLE}"\\t

---> Create list with two columns and use it as input in a do-loop.
while read C ; do
    C1=$(echo ${C} | cut -f1)
    C2=$(echo ${C} | cut -f2)
done <<< "$(paste -d'"'"'\t'"'"' <(echo ${LIST1} | grep -o .) <(echo ${LIST2} | grep -o .))"

---> Insert string at first line of file.
sed -i '"'"'1i<STRING>'"'"' <FILE>

---> Find maximum value in column 3 and print corresponding value in column 2:
<INPUT> | awk -v MAX=0 '"'"'{if($3>MAX){OUTVALUE=$2; MAX=$3}}END{print OUTVALUE}'"'"'

---> Find minimum value in column 3 and print corresponding value of column 2:
<INPUT> | awk -v MIN=0 '"'"'{if($3<MIN){OUTVALUE=$2; MIN=$3}}END{print OUTVALUE}'"'"'

---> Calculate difference between consecutive values in column 2:
<INPUT> | awk '"'"'{if(NR>1){print _n-$2};_n=$2}'"'"'

---> Calculate ratio between two variables:
awk "BEGIN {print ${VAR1}/${VAR2}}"

---> Delete last line.
<INPUT> | sed '"'"'$ d'"'"'

---> Select line that contain a pattern and print the target line and 5 lines before and 5 lines after the target line.
<INPUT> | grep -B 5 -A 2 '"'"'PATTERN'"'"'

---> Find lines where column 2 match an if condition (value below threshold) and print columns 1, 2 and YES/NO (if-else statement with awk).
<INPUT> | awk -v THRESHOLD="${THRESHOLD}" '"'"'{ if ($2 < THRESHOLD) { print $1"\t"$2"\tYES" } else { print $1"\t"$2"\tNO" } }'"'"'

---> Print tab-separated sequence of numbered strings.
printf '"'"'%s\t'"'"' STRING{1..20}

---> Change permissions of files recursively:
find . -type f -exec chmod 644 -- {} +
chmod -R -x+X . #+X gives execute permission only to directories and files that have other execute permissions.

---> Remove all-zeroes rows:
cat ${INPUTFILE} | awk '"'"'{for(i=2;i<=NF;i++) if($i >= 1){print; next}}'"'"'

---> Remove all lines after pattern.
cat ${INPUTFILE} | sed '"'"'1,/<PATTERN>/d'"'"'

---> Remove all lines before pattern.
cat ${INPUTFILE} | sed '"'"'0,/<PATTERN>/d'"'"'

---> Select rows based on file with list of pattern strings.
cat ${INPUTFILE} | grep -wF -f ${PATTERNFILE}

---> Calculate the sum of the numbers in column 3.
cat ${INPUTFILE} | tail -n+2 | cut -f3 | paste -sd+ | bc

---> Substitute/fill empty lines in column 2 with value from column 1.
cat ${INPUTFILE} | awk -F'"'"'\t'"'"' '"'"'$1 && !$2{ $2=$1 }1'"'"' | tr '"'"' '"'"' '"'"'\t'"'"'

---> Print all columns with AWK:
cat ${INPUTFILE} | awk '"'"'{print $0}'"'"'

---> Print all but first column using AWK.
cat ${INPUTFILE} | awk '"'"'{$1=""; print $0}'"'"'

---> Print all but the first two columns using AWK.
cat ${INPUTFILE} | awk '"'"'{$1=$2=""; print $0}'"'"'

---> Get repeated lines and number of occurences.
cat ${INPUTFILE} | sort | uniq -c | awk '$1 > 1 {print $2"\t"$1}'

---> 


---> 


---> 


---> 


---> 


---> 


---> 


---> 


---> 


---> 


---> 


---> 


---> 


---> 


---> 


---> 


---> IMPORTANT PATHS: 

Home directory: ${HOME}
Project home directory: ${PROJHOME}

My scripts directory: ${PATHTOMYSCRIPTSDIR}
My bash profiles: ${PATHTOMYBASHPROFILES}
My scripts: ${PATHTOMYSCRIPTS}
My submitted scripts: ${PATHTOMYSUBMITTEDSCRIPTS}
My SLURM directory: ${PATHTOMYSLURM}

Project working directory: ${PATHTOPROJWORKINGDIR}
Project output directory: ${PATHTOPROJOUTPUT}
Project TMP directory: ${PATHTOPROJTMP}
Trash directory (used in rm-stdout.sh): ${PATHTOPROJTRASH}

Reference genome (Gal6): ${REFGAL6}
Reference genome (Gal7b): ${REFGAL7B}
Reference genome (GRCh38): ${REFHOMOGRCh38}
Chicken consortium data (.bcf): ${CC5KBCF}
' | less

############################################################################
##END OF BASH SCRIPT...
############################################################################