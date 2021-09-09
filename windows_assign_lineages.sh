#!/bin/bash

set -e

#command line inputs
echo "What is the name of the folder containing your sequences?"
read runname
echo "What reference would you like to use? Please type Cosmo_N or Cosmo_WGS"
read reference

mv -vn $runname/*.fasta $runname/$runname".fasta"

mafft --add $runname/$runname".fasta" --reorder inst/extdata/References/$reference/reference_aligned.fasta > $runname/$runname"_withref.fasta"

#Lineage assignment
Rscript Run/windows_run_assignment.R $runname $reference

rm $runname/$runname"_withref.fasta"
