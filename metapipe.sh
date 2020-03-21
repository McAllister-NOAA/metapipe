#!/bin/bash

#NOTES: Require users to supply raw sequence data files in format "samplename_R[12].fastq.gz".
## These names should be identical to names in the sample metadata file.
## Names will be cleaned of illegal characters (only alphanumeric and underline allowed).
## This cleaning step will also be applied to the sample names in the metadata file (so they are the same).
#
#NO spaces or quotes or other weird characters in file paths
#metapipe.sh must be in PATH
#Pro tip: steps that use memory only use 70% of what max you give
#Will need nt and taxonomydmp - run prepscript beforehand
#Can provide own blastn output on whatever database. Must be formatted -outfmt '6 qseqid pident length staxids'

unset parameterfilepath
unset samplemetafilepath
unset readfolderpath
unset outdirectory
unset optionalUserBLASTResult

workingdirectory=`pwd`

unset metapipedir
unset tempprogdir
tempprogdir=`which metapipe.sh`
metapipedir=`echo $tempprogdir | sed -E 's/\/metapipe\.sh$//'`


##########################################################################################
##
##    Get command line options
##
##########################################################################################
pflag=0
sflag=0
rflag=0
oflag=0
blastflag=FALSE

while getopts ":p:s:r:o:b:" opt; do
  case ${opt} in
    p ) pflag=1
        parameterfilepath=$OPTARG #metapipe_config.txt (see README)
      ;;
    s ) sflag=1
        samplemetafilepath=$OPTARG #Sample metadata file (see README)
      ;;
    r ) rflag=1
        readfolderpath=$OPTARG #Folder with fastq.gz raw read files
      ;;
    o ) oflag=1
        outdirectory=$OPTARG #Location for output files
      ;;
    b ) blastflag=TRUE
        optionalUserBLASTResult=$OPTARG #Location of user blastn input (optional)
      ;;
    \? ) echo "Invalid option: -$OPTARG"
         echo "Usage: metapipe.sh" #Invalid option provided
         echo "       -p Config File"
         echo "       -s Sample metadata file"
         echo "       -r Read folder"
         echo "       -o Output directory"
         echo "       -b User-supplied BLASTn btab result file (optional)"
         echo "       See README for details"
         exit
      ;;
    : ) echo "Option is missing an argument: -$OPTARG"
        echo "Usage: metapipe.sh" #Arg for a called option not provided
        echo "       -p Config File"
        echo "       -s Sample metadata file"
        echo "       -r Read folder"
        echo "       -o Output directory"
        echo "       -b User-supplied BLASTn btab result file (optional)"
        echo "       See README for details"
        exit
      ;;
  esac
done
shift $((OPTIND -1))

if [ $OPTIND -eq 1 ]
  then echo "Usage: metapipe.sh" #No options passed
        echo "       -p Config File"
        echo "       -s Sample metadata file"
        echo "       -r Read folder"
        echo "       -o Output directory"
        echo "       -b User-supplied BLASTn btab result file (optional)"
        echo "       See README for details"
        exit
    fi

if [[ $pflag -eq 0 || $sflag -eq 0 || $rflag -eq 0 || $oflag -eq 0 ]]
  then echo "All options except -b are required."
        echo "Usage: metapipe.sh" #Missing required options
        echo "       -p Config File"
        echo "       -s Sample metadata file"
        echo "       -r Read folder"
        echo "       -o Output directory"
        echo "       -b User-supplied BLASTn btab result file (optional)"
        echo "       See README for details"
        exit
      fi
##########################################################################################
##########################################################################################

echo
echo "Start of run:"
date
echo

##########################################################################################
##
##    Get program config file options
##
##########################################################################################
unset primerF
unset primerR
unset ampliconSize
unset desiredTaxaDepth
unset systemmemoryMB
unset dada_minlength
unset dada_phix
unset dada_trunQ
unset dada_maxEE1
unset dada_maxEE2
unset dada_trimRight
unset dada_trimLeft
unset forceMerge
unset deseqRarefaction
unset filterOutContaminants
unset filterInTaxaOfInterest
unset controlPos
unset controlNeg
unset replicates
unset chemData
unset locationNTdatabase
unset speciesGenusCutoffs

source $parameterfilepath

##########################################################################################
##########################################################################################

##########################################################################################
##
##    Cleanup and prep
##
##########################################################################################

mkdir -p ${outdirectory}
mkdir ${outdirectory}/cutadapt
mkdir ${outdirectory}/dada2

#Clean raw read files to remove illegal characters
cd ${readfolderpath}
for f in *fastq.gz
do cleanprior=`echo $f | sed -E 's/^MP_//'`
        base=$(basename $cleanprior .fastq.gz)
        newname=`echo ${base} | sed -E 's/[^A-Za-z0-9_]/_/g'`
        mv $f MP_${newname}.fastq.gz
        done
cd ${workingdirectory}

#Create ordered sample name file
cat ${samplemetafilepath} | cut -f1 | grep -v "Sample" | sed -E 's/[^A-Za-z0-9_]/_/g' | sed -E 's/^/MP_/' > ${outdirectory}/sample_order.txt

##########################################################################################
##########################################################################################

##########################################################################################
##
##    Run CUTADAPT (v.2.8)
##
##########################################################################################
echo "Running Cutadapt: `date`"

#Reverse complement primers (I is converted to N; T/U ok)
revcomp_primerF=`echo $primerF | tr ACGTUWSMKRYBDHVNIacgtuwsmkrybdhvni TGCAAWSKMYRVHDBNNtgcaawskmyrvhdbnn | rev`
revcomp_primerR=`echo $primerR | tr ACGTUWSMKRYBDHVNIacgtuwsmkrybdhvni TGCAAWSKMYRVHDBNNtgcaawskmyrvhdbnn | rev`

#Run cutadapt
for sample in $(cat ${outdirectory}/sample_order.txt)
do
    echo "======================================" >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
    echo "On sample: $sample" >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
    echo "======================================" >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
    echo >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
    cutadapt -a "${primerF};required...${revcomp_primerR};optional" \
    -A "${primerR};required...${revcomp_primerF};optional" \
    --discard-untrimmed \
    -j 0 \
    -m 1 \
    -o ${outdirectory}/cutadapt/${sample}_R1_trimmed.fq.gz -p ${outdirectory}/cutadapt/${sample}_R2_trimmed.fq.gz \
    ${readfolderpath}/${sample}_R1.fastq.gz ${readfolderpath}/${sample}_R2.fastq.gz \
    >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt 2>&1
    echo >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
    echo >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
done

echo "Finished Cutadapt: `date`"

#print cutadapt stats
echo "Sample	Passing Reads	Passing bp"
paste ${outdirectory}/sample_order.txt <(grep "passing" ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")")

echo
echo "Please check Cutadapt success. Proceed? [y/n]"
read mainmenuinput
if [[ "$mainmenuinput" = "y" || "$mainmenuinput" = "Y" ]]; then
  echo "Continuing!"
elif [[ "$mainmenuinput" = "n" || "$mainmenuinput" = "N" ]]; then
  echo "You have chosen to exit"
  exit
else
  echo "Invalid selection, exiting"
  exit
fi

##########################################################################################
##
##    Run DADA2
##
##########################################################################################
echo
echo "Running DADA2: `date`"
echo "Trim and filter in DADA2..."

dada2continue=FALSE
while [ $dada2continue = FALSE ]
  do
  Rscript --vanilla ${metapipedir}/assets/dada2_step1.R ${workingdirectory}/${outdirectory}/dada2 $dada_minlength $dada_phix $dada_trunQ $dada_maxEE1 $dada_maxEE2 $dada_trimRight $dada_trimLeft \
    1>> ${workingdirectory}/${outdirectory}/dada2/dada2_rscripts_out.log 2>&1
  echo
  echo "DADA2 Filtering results:"
  echo "Sample	% Reads Passing"
  awk 'NR>1 { print $1, 100 * ( $3 / $2 ) }' ${workingdirectory}/${outdirectory}/dada2/filtered_out_stats.txt
  echo
  echo "Parameters to modify:"
  echo "minLen,rm.phix,truncQ,maxEE-primer1,maxEE-primer2,trimRight,trimLeft"
  echo "Current settings:"
  echo "${dada_minlength},${dada_phix},${dada_trunQ},${dada_maxEE1},${dada_maxEE2},${dada_trimRight},${dada_trimLeft}"
  echo "Please check DADA2 filtering success. Proceed? [y/n/m]"
  read mainmenuinput
  if [[ "$mainmenuinput" = "y" || "$mainmenuinput" = "Y" ]]; then
    echo "Continuing!"
    dada2continue=TRUE
  elif [[ "$mainmenuinput" = "m" || "$mainmenuinput" = "M" ]]; then
    echo "You have chosen to redo DADA2 filtering with modified settings."
    echo "Input new settings separated by commas in the same order as above."
    read secondmainmenuinput
    IFS=','
    read -ra ADDR <<< "$secondmainmenuinput"
    dada_minlength=${ADDR[0]}
    dada_phix=${ADDR[1]}
    dada_trunQ=${ADDR[2]}
    dada_maxEE1=${ADDR[3]}
    dada_maxEE2=${ADDR[4]}
    dada_trimRight=${ADDR[5]}
    dada_trimLeft=${ADDR[6]}
    echo "New settings:"
    echo "${dada_minlength},${dada_phix},${dada_trunQ},${dada_maxEE1},${dada_maxEE2},${dada_trimRight},${dada_trimLeft}"
    echo "Rerunning DADA2 filtering"
  elif [[ "$mainmenuinput" = "n" || "$mainmenuinput" = "N" ]]; then
    echo "You have chosen to exit"
    exit
  else
    echo "Invalid selection, try again"
  fi
done

echo
echo "Learning error, Dereplication, Merge, and ASVs in DADA2..."
echo "Please be patient, may take a while. Messages printed to Rscript log."
echo
Rscript --vanilla ${metapipedir}/assets/dada2_step2.R ${workingdirectory}/${outdirectory}/dada2 $systemmemoryMB \
    1>> ${workingdirectory}/${outdirectory}/dada2/dada2_rscripts_out.log 2>&1

cat ${outdirectory}/dada2/ASVs_counts.tsv | sed -E 's/^	/x	/' > ${outdirectory}/dada2/ASVs_counts_mod.tsv
mv ${outdirectory}/dada2/ASVs_counts_mod.tsv ${outdirectory}/dada2/ASVs_counts.tsv

##TODO:
#HOW TO FORCE CONTINUE WITHOUT MERGE?
#should we allow option for changing minOverlap? (right now set to 20bp, which seems reasonable, don't want chimeras)
#Optional rerun with deseq for rarefaction

echo
echo "FINAL DADA2 STATS"
echo "Sample	%Reads Retained"
awk 'NR>1 { print $1, $8 }' ${workingdirectory}/${outdirectory}/dada2/ReadTrimSummary.txt
echo
echo "Do you wish to Proceed? [y/n]"
read mainmenuinput
if [[ "$mainmenuinput" = "y" || "$mainmenuinput" = "Y" ]]; then
  echo "Continuing!"
elif [[ "$mainmenuinput" = "n" || "$mainmenuinput" = "N" ]]; then
  echo "You have chosen to exit"
  exit
else
  echo "Invalid selection, exiting"
  exit
fi

##########################################################################################
##
##    Blastn ASVs
##
##########################################################################################
mkdir ${outdirectory}/blast_results

if [[ "$blastflag" = "FALSE" ]]; then
echo
echo "Running BLASTn: `date`"

passblastscrutiny=FALSE
maxtargetseqs=500
runthroughcount=0

while [ $passblastscrutiny = FALSE ]
  do
  let "runthroughcount=runthroughcount+1"
  blastn -db ${locationNTdatabase}/nt -query ${outdirectory}/dada2/ASVs.fa -outfmt '6 qseqid pident length staxids' -max_target_seqs $maxtargetseqs -num_threads 6 -out ${outdirectory}/blast_results/ASV_blastn_nt.btab
  blastn -db ${locationNTdatabase}/nt -query ${outdirectory}/dada2/ASVs.fa -outfmt '6 qseqid pident length staxids' -max_target_seqs $maxtargetseqs -num_threads 6 -out ${outdirectory}/blast_results/ASV_blastn_nt_ignoreMostEnvSeqs_leavesinUnclassified.btab -negative_taxidlist ${locationNTdatabase}/taxdump/taxid_exclusion_list_leavesinUnclassified.txt
  blastn -db ${locationNTdatabase}/nt -query ${outdirectory}/dada2/ASVs.fa -outfmt '6 qseqid pident length staxids' -max_target_seqs $maxtargetseqs -num_threads 6 -out ${outdirectory}/blast_results/ASV_blastn_nt_ignoreAllEnvSeqs_removesUnclassified.btab -negative_taxidlist ${locationNTdatabase}/taxdump/taxid_exclusion_list_removesUnclassified.txt
  for f in $(eval echo "{1..`grep -c ">" ${outdirectory}/dada2/ASVs.fa`}")
    do echo -ne "ASV_${f}\t" >> ${outdirectory}/blast_results/fullhits_checkmaxtargetseqs.txt
      if grep -q "ASV_${f}\t" ${outdirectory}/blast_results/ASV_blastn_nt.btab; then
        grep "ASV_${f}\t" ${outdirectory}/blast_results/ASV_blastn_nt.btab | cut -f2 | sort -nr | uniq -c | head -1 | sed -E 's/^ +//' | sed -E 's/([0-9]) ([0-9])/\1	\2/' >> ${outdirectory}/blast_results/fullhits_checkmaxtargetseqs.txt
      else echo "0 NONE" >> ${outdirectory}/blast_results/fullhits_checkmaxtargetseqs.txt
      fi
  done
  for f in $(eval echo "{1..`grep -c ">" ${outdirectory}/dada2/ASVs.fa`}")
    do echo -ne "ASV_${f}\t" >> ${outdirectory}/blast_results/ignoreMostEnvSeqs_checkmaxtargetseqs.txt
      if grep -q "ASV_${f}\t" ${outdirectory}/blast_results/ASV_blastn_nt_ignoreMostEnvSeqs_leavesinUnclassified.btab; then
        grep "ASV_${f}\t" ${outdirectory}/blast_results/ASV_blastn_nt_ignoreMostEnvSeqs_leavesinUnclassified.btab | cut -f2 | sort -nr | uniq -c | head -1 | sed -E 's/^ +//' | sed -E 's/([0-9]) ([0-9])/\1	\2/' >> ${outdirectory}/blast_results/ignoreMostEnvSeqs_checkmaxtargetseqs.txt
      else echo "0 NONE" >> ${outdirectory}/blast_results/ignoreMostEnvSeqs_checkmaxtargetseqs.txt
      fi
  done
  for f in $(eval echo "{1..`grep -c ">" ${outdirectory}/dada2/ASVs.fa`}")
    do echo -ne "ASV_${f}\t" >> ${outdirectory}/blast_results/ignoreAllEnvSeqs_checkmaxtargetseqs.txt
      if grep -q "ASV_${f}\t" ${outdirectory}/blast_results/ASV_blastn_nt_ignoreAllEnvSeqs_removesUnclassified.btab; then
        grep "ASV_${f}\t" ${outdirectory}/blast_results/ASV_blastn_nt_ignoreAllEnvSeqs_removesUnclassified.btab | cut -f2 | sort -nr | uniq -c | head -1 | sed -E 's/^ +//' | sed -E 's/([0-9]) ([0-9])/\1	\2/' >> ${outdirectory}/blast_results/ignoreAllEnvSeqs_checkmaxtargetseqs.txt
      else echo "0 NONE" >> ${outdirectory}/blast_results/ignoreAllEnvSeqs_checkmaxtargetseqs.txt
      fi
  done
  for f in ${outdirectory}/blast_results/*checkmaxtargetseqs.txt
    do cat $f | sort -k2 -nr | head -1 | cut -f2 >> ${outdirectory}/blast_results/temp.txt
  done
  highestnumber=`cat ${outdirectory}/blast_results/temp.txt | sort -nr | head -1`
  rm ${outdirectory}/blast_results/temp.txt
  if [[ "$highestnumber" -lt "$maxtargetseqs" ]]; then
    passblastscrutiny=TRUE
  else
    echo "Rerun BLAST Number ${runthroughcount} - Max Target Not High Enough (${maxtargetseqs})"
    let "maxtargetseqs=maxtargetseqs*2"
    echo "New max target seqs = ${maxtargetseqs}"
    rm ${outdirectory}/blast_results/fullhits_checkmaxtargetseqs.txt
    rm ${outdirectory}/blast_results/ignoreMostEnvSeqs_checkmaxtargetseqs.txt
    rm ${outdirectory}/blast_results/ignoreAllEnvSeqs_checkmaxtargetseqs.txt
  fi
done

elif [[ "$blastflag" = "TRUE" ]]; then
cp $optionalUserBLASTResult ${outdirectory}/blast_results/ASV_blastn_USERSUPPLIED.btab
fi

##########################################################################################
##
##    Reformat BLAST output
##
##########################################################################################
echo
echo "Reformatting BLAST output: `date`"

let "insertSize=ampliconSize - ${#primerF} - ${#primerR}"

if [[ "$blastflag" = "FALSE" ]]; then
Rscript --vanilla ${metapipedir}/assets/reformat_blast_default.R ${workingdirectory}/${outdirectory}/blast_results $insertSize \
    1>> ${workingdirectory}/${outdirectory}/blast_results/blastreformatting_rscript_out.log 2>&1
elif [[ "$blastflag" = "TRUE" ]]; then
Rscript --vanilla ${metapipedir}/assets/reformat_blast_user.R ${workingdirectory}/${outdirectory}/blast_results $insertSize \
    1>> ${workingdirectory}/${outdirectory}/blast_results/blastreformatting_userprovided_rscript_out.log 2>&1
echo
fi

##########################################################################################
##
##    Running perl script for taxonomy decisions, standard table outs, and KRONA plots
##
##########################################################################################
echo
echo "Running ASV-2-Taxonomy Script: `date`"

mkdir ${outdirectory}/ASV2Taxonomy

#prep taxon ids for taxonkit
cat ${outdirectory}/blast_results/*formatted.txt | cut -f4 | sed -E "s/ //g" | tr ',' '\n' | tr ';' '\n' | sort | uniq | grep -v "taxid" > ${outdirectory}/ASV2Taxonomy/taxids.txt
taxonkit lineage ${outdirectory}/ASV2Taxonomy/taxids.txt | awk '$2!=""' > ${outdirectory}/ASV2Taxonomy/taxonkit_out.txt
taxonkit reformat ${outdirectory}/ASV2Taxonomy/taxonkit_out.txt | cut -f1,3 > ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt

echo
echo "Reformatted taxon strings created. Check and edit if you wish. Proceed? [any]"
read mainmenuinput
if [[ "$mainmenuinput" = "y" ]]; then
  echo "Continuing!"
else
  echo "Continuing!"
fi


${metapipedir}/assets/convert_line_endings.py ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt #Deprecated 'U' option. TODO find replacement

cd ${outdirectory}/ASV2Taxonomy
perl ${metapipedir}/assets/asv_taxonomy_processing_figureOuts.pl -a ../dada2/ASVs_counts.tsv -s ../blast_results/ASV_blastn_nt_formatted.txt \
    -t reformatted_taxonkit_out.txt -f $speciesGenusCutoffs -n ${outdirectory} -c ${locationNTdatabase}/taxdump/common_names.dmp -o ../sample_order.txt

cd ${workingdirectory}


#TODO will need to decide on what blastn outfiles to focus on and what to use. For now, just running all here... Need to discuss.












echo "YOU MADE IT!"









