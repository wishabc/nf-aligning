#! /bin/bash

HOTSPOT_DISTR=$1
OUTDIR=$2
TAGS=$3
GENOME=$4
chromInfo=$5
K=$6
ASSAY=$7

echo "Making SPOT in $OUTDIR"
echo "TAGS: $TAGS"
echo "GENOME: $GENOME"
echo "K: $K"
echo "ASSAY: $ASSAY"

CONFIGOUT=$OUTDIR/runall.tokens.txt
SCRIPTOUT=$OUTDIR/runhotspot

if [ ! -e $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

shopt -s nocasematch
# Decide whether to use DUPOK based on assay type
if [[ "$ASSAY" == *"ChIP"* ]]
then
    echo "Setting DUPOK=F"
    DUPOK="F"
else
    echo "Setting DUPOK=T"
    DUPOK="T"
fi

echo "Creating $CONFIGOUT"

chromFile=$(readlink -f "$5")
mappableFile=$(readlink -f "$GENOME.K$K.mappable_only.bed")
chkchr=$(awk 'NR==1{print $1}' < "$chromFile")

# Create the configuration file for the SPOT program
cat > $CONFIGOUT <<EOF
[script-tokenizer]

#######################################
## Notes:  If duplicate token definitions exist, the last definition in
##         the file will be used. Tokens can use .ini variables and declarations.
##         See http://docs.python.org/library/configparser.html
#######################################

#######################################
# Global tokens (used by most scripts)
#######################################

## Tags file in bam format (file extension .bam), or starched bed file
## (file extension .bed.starch).  If the tags file is in bed.starch
## format, and if it is in directory specified by _OUTDIR_, then it
## needs to be in the format
##
##     chr  5'start  5'start+1
##
## That is, the file should be three column bed (no ID field, etc.)
## containing the 1bp coordinates of the 5' ends of the tags. If
## _TAGS_ is a bam file, or a bed.starch file not in _OUTDIR_, then
## the bed.starch file in the above format will be generated (using
## the strand column if present), and put in _OUTDIR_.  NOTE: if you
## use run_badspot, then you must use bam files, or a bed.starch file
## not in _OUTDIR_.

_TAGS_ = $TAGS

## For ChIP data with an Input bam file, set _USE_INPUT_ to T, and set
## _INPUT_TAGS_ to the name of that bam file.

_USE_INPUT_ = F

## Genome
_GENOME_ = $GENOME
## Tag length
_K_ = $K

## Chromosome coordinates, bed format.
_CHROM_FILE_ = $chromFile
## Location of uniquely mappable positions in the genome for this tag length.
_MAPPABLE_FILE_ = $mappableFile

## Set DUPOK to T for DNaseI data, F for ChIP-seq data (DUPOK = T means allow duplicate reads)
_DUPOK_ = $DUPOK

## FDR levels, separated by spaces if more than one. Set to N if you
## do not want FDR thresholding (for example, if you just want SPOT
## score computed.)
_FDRS_ = "0.05,0.01,0.001"

## Output directories (can all be the same location).  Use full path names.
## _OUTDIR_ contains tags files in converted bed.starch and lib.txt formats (for hotspot
## program), and hotspot and peak results.
## _RANDIR_ contains generated random tags (for FDR thresholding) and hotspots called on random tags.
_OUTDIR_ = $OUTDIR
_RANDIR_ = $OUTDIR

## Set to T if you want scripts to skip steps that have already been done.
_CHECK_ = T

## If _CHECK_ = T, outputs are checked for completeness by searching
## for results for the following chromsome.
_CHKCHR_ = $chkchr

## Hotspot program binary
_HOTSPOT_ = $HOTSPOT_DISTR/hotspot-deploy/bin/hotspot

## Clean up. Remove all intermediate files and directories if set to T.  See
## pipeline script run_final.
_CLEAN_ = T

_PKFIND_BIN_ = $HOTSPOT_DISTR/hotspot-deploy/bin/wavePeaks

## Peak-finding smoothing level. If the resolution of the input file
## is x, then the results are smoothed out to a scale of (2^level)*x.
_PKFIND_SMTH_LVL_ = 3

## Random number seed, used for generating random tags for FDR thresholding.
_SEED_=101

## Hotspot program parameters
_THRESH_ = 2
_WIN_MIN_ = 200
_WIN_MAX_ = 300
_WIN_INCR_ = 50
_BACKGRD_WIN_ = 50000
_MERGE_DIST_ = 150
_MINSIZE_ = 10
EOF

# Save the basic script we are running for easy rerunning 
# in troubleshooting
cat > $SCRIPTOUT <<EOF

scriptTokBin="python2 ${HOTSPOT_DISTR}/ScriptTokenizer/src/script-tokenizer.py"
pipeDir=${HOTSPOT_DISTR}/pipeline-scripts
tokenFile=$CONFIGOUT

cd $OUTDIR

# Do SPOT only (set _FDRS_ to "N" in runall.tokens.txt)

scripts="
    \$pipeDir/run_badspot
    \$pipeDir/run_make_lib
    \$pipeDir/run_wavelet_peak_finding
    \$pipeDir/run_10kb_counts
    \$pipeDir/run_generate_random_lib
    \$pipeDir/run_pass1_hotspot
    \$pipeDir/run_pass1_merge_and_thresh_hotspots
    \$pipeDir/run_pass2_hotspot
    \$pipeDir/run_rescore_hotspot_passes
    \$pipeDir/run_spot
    \$pipeDir/run_thresh_hot.R
    \$pipeDir/run_both-passes_merge_and_thresh_hotspots
    \$pipeDir/run_add_peaks_per_hotspot
    \$pipeDir/run_final
"

\$scriptTokBin \
    --clobber \
    --output-dir=$OUTDIR \
    \$tokenFile \
    \$scripts

wavelets=$HOTSPOT_DISTR/hotspot-deploy/bin/wavelets
which wavelets
for script in \$scripts
do
    ./\$(basename \$script).tok
done
EOF

echo "Running $SCRIPTOUT"

bash $SCRIPTOUT