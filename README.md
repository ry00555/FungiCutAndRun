# FungiCutAndRun

#number of CPUs
THREADS=12

#FASTQ should be the path to a folder containing all your zipped fastq files; filename endings should be fastq.gz

FASTQ="path/to/fastqfiles"

#OUTDIR should be the path a folder for writing output
OUTDIR="path/to/OutputDirectory"

#GENOME should be the path to a bwa-indexed genome file; include the entire root file name
GENOME="path/to/IndexedGenomeFile"

#SET Deeptools paramenters
#binsize for windows
BIN="25"

#smoothlength setting, for smoothing the enrichment curve
$SMOOTH="50"

#if you want to use the MNase parameter, uncomment the following lines
MNase="--MNase"
MN="MNase"

