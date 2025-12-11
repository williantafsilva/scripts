#!/bin/bash -l
############################################################################
############################### TEST SCRIPT ################################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

module load PDC/24.11
module load R/4.4.2-cpeGNU-24.11
INPUTFILE_CIS=$(readlink -f matrixeqtl-ciseqtl-noFDR-job14602496.txt)
INPUTFILE_TRANS=$(readlink -f matrixeqtl-transeqtl-noFDR-job14602496.txt)
PVALCOLUMN=5
NTESTS_CIS=20335031958
NTESTS_TRANS=10637201907434
paste ${INPUTFILE_CIS} <(Pvalue-FDR-stdout.sh ${INPUTFILE_CIS} ${PVALCOLUMN} ${NTESTS_CIS} | sed '1iFDR') > matrixeqtl-ciseqtl-FDR.txt
paste ${INPUTFILE_TRANS} <(Pvalue-FDR-stdout.sh ${INPUTFILE_TRANS} ${PVALCOLUMN} ${NTESTS_TRANS} | sed '1iFDR') > matrixeqtl-transeqtl-FDR.txt

touch FDRcalculationfinished