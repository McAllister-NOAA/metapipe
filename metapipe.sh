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

#TODO will need a better catch checking if sample names from metadata file match read sample names.


unset parameterfilepath
unset samplemetafilepath
unset readfolderpath
unset outdirectory
unset optionalUserBLASTResult
unset figureparamfilepath

workingdirectory=`pwd`

unset metapipedir
unset tempprogdir
tempprogdir=`which metapipe.sh`
metapipedir=`echo $tempprogdir | sed -E 's/\/metapipe\.sh$//'`
myInvocation="$(printf %q "$BASH_SOURCE")$((($#)) && printf ' %q' "$@")"

##########################################################################################
##
##    Get command line options
##
##########################################################################################
pflag=0
sflag=0
rflag=0
oflag=0
fflag=0
blastflag=FALSE
bypassflag=FALSE

while getopts ":p:s:r:o:b:f:y" opt; do
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
        outdirectory=`echo $outdirectory | sed -E 's/\/$//'`
      ;;
    f ) fflag=1
        figureparamfilepath=$OPTARG #Location of figure config file (see README)
      ;;
    b ) blastflag=TRUE
        optionalUserBLASTResult=$OPTARG #Location of user blastn input (optional)
      ;;
    y ) bypassflag=TRUE
      ;;
    \? ) echo "Invalid option: -$OPTARG"
         echo "Usage: metapipe.sh" #Invalid option provided
         echo "       -p Config File"
         echo "       -f Figure config file"
         echo "       -s Sample metadata file"
         echo "       -r Read folder"
         echo "       -o Output directory"
         echo "       -b User-supplied BLASTn btab result file (optional)"
         echo "       -y Bypass all terminal prompts (optional)"
         echo "       See README for details"
         exit
      ;;
    : ) echo "Option is missing an argument: -$OPTARG"
        echo "Usage: metapipe.sh" #Arg for a called option not provided
        echo "       -p Config File"
        echo "       -f Figure config file"
        echo "       -s Sample metadata file"
        echo "       -r Read folder"
        echo "       -o Output directory"
        echo "       -b User-supplied BLASTn btab result file (optional)"
        echo "       -y Bypass all terminal prompts (optional)"
        echo "       See README for details"
        exit
      ;;
  esac
done
shift $((OPTIND -1))

if [ $OPTIND -eq 1 ]
  then echo "Usage: metapipe.sh" #No options passed
        echo "       -p Config File"
        echo "       -f Figure config file"
        echo "       -s Sample metadata file"
        echo "       -r Read folder"
        echo "       -o Output directory"
        echo "       -b User-supplied BLASTn btab result file (optional)"
        echo "       -y Bypass all terminal prompts (optional)"
        echo "       See README for details"
        exit
    fi

if [[ $pflag -eq 0 || $sflag -eq 0 || $rflag -eq 0 || $oflag -eq 0 || $fflag -eq 0 ]]
  then echo "All options except -b and -y are required."
        echo "Usage: metapipe.sh" #Missing required options
        echo "       -p Config File"
        echo "       -f Figure config file"
        echo "       -s Sample metadata file"
        echo "       -r Read folder"
        echo "       -o Output directory"
        echo "       -b User-supplied BLASTn btab result file (optional)"
        echo "       -y Bypass all terminal prompts (optional)"
        echo "       See README for details"
        exit
      fi
##########################################################################################
##########################################################################################

##########################################################################################
##
##    Get program config file options
##
##########################################################################################
unset primerF
unset primerR
unset ampliconSize
unset systemmemoryMB
unset locationNTdatabase
unset speciesGenusCutoffs
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
unset blastMode

source $parameterfilepath

##########################################################################################
##########################################################################################

##########################################################################################
##
##    Cleanup and prep
##
##########################################################################################
cutadaptFinished=FALSE
dada2_Finished=FALSE
blastFinished=FALSE
taxonomyscriptFinished=FALSE
figuresFinished=FALSE

if [ -d "${outdirectory}" ]; then
  configcompare1=`cat ${parameterfilepath}`
  configcompare2=`cat ${outdirectory}/config_file.txt`
  samplecompare1=`cat ${samplemetafilepath}`
  samplecompare2=`cat ${outdirectory}/sample_metadata.txt`
  if [[ "${configcompare1}" = "${configcompare2}" ]]; then
    echo "Config files identical"
  else
    echo "Config files differ between runs, choose a different out directory"
    exit
  fi
  if [[ "${samplecompare1}" = "${samplecompare2}" ]]; then
    echo "Sample metadata files identical"
  else
    echo "Sample metadata files differ between runs, choose a different out directory"
    exit
  fi
  source ${outdirectory}/progress.txt
  mkdir -p ${outdirectory}/run_logs
  currenttime=`date | sed -E 's/[^A-Za-z0-9_]/_/g'`
  mv ${outdirectory}/run.log ${outdirectory}/run_logs/runlog_${currenttime}.txt
  touch ${outdirectory}/run.log
else
  mkdir ${outdirectory}
  touch ${outdirectory}/progress.txt
  touch ${outdirectory}/run.log
  cp ${parameterfilepath} ${outdirectory}/config_file.txt
  cp ${samplemetafilepath} ${outdirectory}/sample_metadata.txt
fi

exec &> >(tee -a ${outdirectory}/run.log)

echo
echo "Start of run:"
date
echo
echo "Invoked script options:"
echo "$myInvocation"
echo

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

#Create sample metadata file with identical manipulation of sample names for downstream R work
cat ${samplemetafilepath} | awk 'FNR==NR {gsub("[^a-zA-Z0-9]", "_", $1)} 1' OFS="\t" | sed -e '2,$ s/^/MP_/' | grep -v "^MP_$" > ${outdirectory}/sample_metadata_forR.txt

##########################################################################################
##########################################################################################

##########################################################################################
##
##    Run CUTADAPT (v.2.8)
##
##########################################################################################
if [[ "${cutadaptFinished}" = "TRUE" ]]; then
  echo "Cutadapt from prior run"
else
  echo "Running Cutadapt: `date`"
  if [ -d "${outdirectory}/cutadapt" ]; then
    rm -r ${outdirectory}/cutadapt
    mkdir ${outdirectory}/cutadapt
  else
    mkdir ${outdirectory}/cutadapt
  fi
  
  #Reverse complement primers (I is converted to N; T/U ok)
  revcomp_primerF=`echo $primerF | tr ACGTUWSMKRYBDHVNIacgtuwsmkrybdhvni TGCAAWSKMYRVHDBNNtgcaawskmyrvhdbnn | rev`
  revcomp_primerR=`echo $primerR | tr ACGTUWSMKRYBDHVNIacgtuwsmkrybdhvni TGCAAWSKMYRVHDBNNtgcaawskmyrvhdbnn | rev`
  
  temp_primerF=`echo $primerF | tr Ii Nn`
  temp_primerR=`echo $primerR | tr Ii Nn`
  primerF=$temp_primerF
  primerR=$temp_primerR  
  
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
  
  if [[ "$bypassflag" = "FALSE" ]]; then
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
  fi

  echo "cutadaptFinished=TRUE" >> ${outdirectory}/progress.txt

fi

#TODO: What to do with samples where ALL reads are filtered out? Will it produce an error?
#TODO: May consider removing intermediate fastq.gz files after they are used. They take up a lot of space.

##########################################################################################
##
##    Run DADA2
##
##########################################################################################
if [[ "${dada2_Finished}" = "TRUE" ]]; then
  echo "DADA2 from prior run"
else
  echo
  echo "Running DADA2: `date`"
  if [ -d "${outdirectory}/dada2" ]; then
    rm -r ${outdirectory}/dada2
    mkdir ${outdirectory}/dada2
  else
    mkdir ${outdirectory}/dada2
  fi
  
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
    if [[ "$bypassflag" = "FALSE" ]]; then
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
    else
      dada2continue=TRUE
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
  if [[ "$bypassflag" = "FALSE" ]]; then
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
  fi

  echo "dada2_Finished=TRUE" >> ${outdirectory}/progress.txt

fi

##########################################################################################
##
##    Blastn ASVs
##
##########################################################################################
if [[ "${blastFinished}" = "TRUE" ]]; then
  echo "BLASTn from prior run"
else
  if [ -d "${outdirectory}/blast_results" ]; then
    rm -r ${outdirectory}/blast_results
    mkdir ${outdirectory}/blast_results
  else
    mkdir ${outdirectory}/blast_results
  fi
  
  if [[ "$blastflag" = "FALSE" ]]; then
  echo
  echo "Running BLASTn: `date`"
  
  passblastscrutiny=FALSE
  maxtargetseqs=400
  runthroughcount=0
  
  while [ $passblastscrutiny = FALSE ]
    do
    let "runthroughcount=runthroughcount+1"
    if [[ "${blastMode}" = "allIN" ]]; then
      blastn -db ${locationNTdatabase}/nt -query ${outdirectory}/dada2/ASVs.fa -outfmt '6 qseqid pident length staxids' -subject_besthit -max_target_seqs $maxtargetseqs -num_threads 6 -out ${outdirectory}/blast_results/ASV_blastn_nt.btab
    elif [[ "${blastMode}" = "mostEnvOUT" ]]; then
      blastn -db ${locationNTdatabase}/nt -query ${outdirectory}/dada2/ASVs.fa -outfmt '6 qseqid pident length staxids' -subject_besthit -max_target_seqs $maxtargetseqs -num_threads 6 -out ${outdirectory}/blast_results/ASV_blastn_nt.btab -negative_taxidlist ${locationNTdatabase}/taxdump/taxid_exclusion_list_leavesinUnclassified.txt
    elif [[ "${blastMode}" = "allEnvOUT" ]]; then
      blastn -db ${locationNTdatabase}/nt -query ${outdirectory}/dada2/ASVs.fa -outfmt '6 qseqid pident length staxids' -subject_besthit -max_target_seqs $maxtargetseqs -num_threads 6 -out ${outdirectory}/blast_results/ASV_blastn_nt.btab -negative_taxidlist ${locationNTdatabase}/taxdump/taxid_exclusion_list_removesUnclassified.txt
    else
      echo "Incorrect blastMode specified"
      exit
    fi
    perl ${metapipedir}/assets/blast_assessment.pl -i ${outdirectory}/blast_results/ASV_blastn_nt.btab -c `grep -c ">" ${outdirectory}/dada2/ASVs.fa` > ${outdirectory}/blast_results/checkmaxtargetseqs.txt
    highestnumber=`cat ${outdirectory}/blast_results/checkmaxtargetseqs.txt | sort -k2 -nr | head -1 | cut -f2`
    if [[ "$highestnumber" -lt "$maxtargetseqs" ]]; then
      passblastscrutiny=TRUE
    elif [[ "$highestnumber" -ge "$maxtargetseqs" && "$runthroughcount" -le 3 ]]; then
      echo "Rerun BLAST Number ${runthroughcount} - Max Target Not High Enough (${maxtargetseqs})"
      let "multiplierseqs=runthroughcount+2"
      let "maxtargetseqs=highestnumber*multiplierseqs"
      echo "New max target seqs = ${maxtargetseqs}"
    else
      echo "Max target seqs inclusive of all top hits cannot be reached"
      echo "Used max target seqs parameter = ${maxtargetseqs}"
      passblastscrutiny=TRUE
    fi
  done
  
  elif [[ "$blastflag" = "TRUE" ]]; then
  echo "Using user supplied BLASTn btab result"
  cp $optionalUserBLASTResult ${outdirectory}/blast_results/ASV_blastn_nt.btab
  fi


##########################################################################################
##
##    Reformat BLAST output
##
##########################################################################################
  echo
  echo "Reformatting BLAST output: `date`"
  
  let "insertSize=ampliconSize - ${#primerF} - ${#primerR}"
  
  Rscript --vanilla ${metapipedir}/assets/reformat_blast.R ${workingdirectory}/${outdirectory}/blast_results $insertSize \
      1>> ${workingdirectory}/${outdirectory}/blast_results/blastreformatting_rscript_out.log 2>&1
  
  echo
  
  echo "blastFinished=TRUE" >> ${outdirectory}/progress.txt

fi


##########################################################################################
##
##    Running perl script for taxonomy decisions, standard table outs, and KRONA plots
##
##########################################################################################
if [[ "${taxonomyscriptFinished}" = "TRUE" ]]; then
  echo "ASV-2-Taxonomy Script results from prior run"
else
  if [ -d "${outdirectory}/ASV2Taxonomy" ]; then
    rm -r ${outdirectory}/ASV2Taxonomy
    mkdir ${outdirectory}/ASV2Taxonomy
  else
    mkdir ${outdirectory}/ASV2Taxonomy
  fi

  echo
  echo "Running ASV-2-Taxonomy Script: `date`"
  
  #prep taxon ids for taxonkit
  cat ${outdirectory}/blast_results/ASV_blastn_nt_formatted.txt | cut -f4 | sed -E "s/ //g" | tr ',' '\n' | tr ';' '\n' | sort | uniq | grep -v "taxid" > ${outdirectory}/ASV2Taxonomy/taxids.txt
  taxonkit lineage ${outdirectory}/ASV2Taxonomy/taxids.txt | awk '$2!=""' > ${outdirectory}/ASV2Taxonomy/taxonkit_out.txt
  taxonkit reformat ${outdirectory}/ASV2Taxonomy/taxonkit_out.txt | cut -f1,3 > ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt
  
  if [[ "$bypassflag" = "FALSE" ]]; then
    echo
    echo "Reformatted taxon strings created. Options:"
    echo "Continue without changes [c]"
    echo "Manually edit file and replace in same location with identical file structure [m]"
    echo "    (Make choice when file is modified and you are ready to proceed)"
    echo "Automatically fill gaps in reformatted taxonkit hierarchy [a]"
    read mainmenuinput
    if [[ "$mainmenuinput" = "c" || "$mainmenuinput" = "C" ]]; then
      echo "Continuing!"
    elif [[ "$mainmenuinput" = "m" || "$mainmenuinput" = "M" ]]; then
      echo "Continuing!"
      cat ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt | tr '\r' '\n' | tr -s '\n' > ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp
      mv ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt
    elif [[ "$mainmenuinput" = "a" || "$mainmenuinput" = "A" ]]; then
      echo "Reformatting..."
      echo "Original reformatted taxonkit out stored at ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out_ORIGINAL.txt"
      perl ${metapipedir}/assets/fillIn_taxonkit.pl -i ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt > ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp
      mv ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out_ORIGINAL.txt
      mv ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt
      echo "Continuing!"
    else
      echo "Invalid selection, exiting"
      exit
    fi
  elif [[ "$bypassflag" = "TRUE" ]]; then
    echo "Automatically reformatting taxonkit hierarchy..."
    echo "Original reformatted taxonkit out stored at ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out_ORIGINAL.txt"
    perl ${metapipedir}/assets/fillIn_taxonkit.pl -i ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt > ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp
    mv ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out_ORIGINAL.txt
    mv ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt
    echo "Continuing!"
  fi
    
  cd ${outdirectory}/ASV2Taxonomy
  perl ${metapipedir}/assets/asv_taxonomy_processing_figureOuts.pl -a ../dada2/ASVs_counts.tsv -s ../blast_results/ASV_blastn_nt_formatted.txt \
      -t reformatted_taxonkit_out.txt -f $speciesGenusCutoffs -n ${outdirectory} -c ${locationNTdatabase}/taxdump/common_names.dmp -o ../sample_order.txt
  cat ${outdirectory}_asvTaxonomyTable.txt | grep -v "Unknown" > ${outdirectory}_asvTaxonomyTable_NOUNKNOWNS.txt
  #TODO: cleanup outputs and directory from my asv to taxonomy script
  #TODO: it would be nice to get some stats on how many ASVs could be classified to species, genus, family, etc.
  cd ${workingdirectory}
  cat ${outdirectory}/ASV2Taxonomy/${outdirectory}_unknown_asvids.txt | cut -f1 | sed -E 's/$/	/' > ${outdirectory}/ASV2Taxonomy/temp_grep_unknowns
  cat ${outdirectory}/dada2/ASVs_counts.tsv | grep -v -f ${outdirectory}/ASV2Taxonomy/temp_grep_unknowns > ${outdirectory}/ASV2Taxonomy/ASVs_counts_NOUNKNOWNS.tsv
  perl ${metapipedir}/assets/merge_on_taxonomy.pl -a ${outdirectory}/dada2/ASVs_counts.tsv -t ${outdirectory}/ASV2Taxonomy/${outdirectory}_asvTaxonomyTable.txt > ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy.tsv
  cat ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy.tsv | grep -v -f ${outdirectory}/ASV2Taxonomy/temp_grep_unknowns > ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy_NOUNKNOWNS.tsv
  rm ${outdirectory}/ASV2Taxonomy/temp_grep_unknowns
  mkdir ${outdirectory}/ASV2Taxonomy/KRONA_plots
  mkdir ${outdirectory}/ASV2Taxonomy/KRONA_plots/KRONA_inputs
  mv ${outdirectory}/ASV2Taxonomy/MP* ${outdirectory}/ASV2Taxonomy/KRONA_plots/KRONA_inputs/
  mv ${outdirectory}/ASV2Taxonomy/${outdirectory}_wholeKRONA.txt ${outdirectory}/ASV2Taxonomy/KRONA_plots/KRONA_inputs/${outdirectory}_samplesSummedKRONA.txt
  mv ${outdirectory}/ASV2Taxonomy/${outdirectory}_master_krona.html ${outdirectory}/ASV2Taxonomy/KRONA_plots/
  mv ${outdirectory}/ASV2Taxonomy/${outdirectory}_wholeKRONA.html ${outdirectory}/ASV2Taxonomy/KRONA_plots/${outdirectory}_samplesSummedKRONA.html
  mv ${outdirectory}/ASV2Taxonomy/${outdirectory}_allin_KRONA.txt ${outdirectory}/ASV2Taxonomy/${outdirectory}_allin_TaxonomyASVSampleCount_byline.txt
  
  echo "taxonomyscriptFinished=TRUE" >> ${outdirectory}/progress.txt

fi

##########################################################################################
##
##    File cleanup
##
##########################################################################################
rm -f ${outdirectory}/cutadapt/*_trimmed.fq.gz
gzip -q -9 ${outdirectory}/blast_results/ASV_blastn_nt.btab

##########################################################################################
##
##    FIGURES & ANALYSIS
##
##########################################################################################
#Defaults, user can override:
filterPercent=5
filterLowQualSamples="FALSE"
filterPercentLowQualSamples=30
removeNA="FALSE"
providedTaxaOfInterest="FALSE"
taxaOfInterestFile="NULL"
taxaOfInterestLevel="NULL"
#Interpreted from sample metadata file:
controlPos="FALSE"
controlNeg="FALSE"
replicates="FALSE"
sites="FALSE"
chemData="FALSE"
controlspresent="FALSE"
groupsDefinedFlag="FALSE"
locationChemHeaders="NULL"

source $figureparamfilepath

cp $figureparamfilepath ${outdirectory}/figure_config_file.txt

#TODO Check previous and current config files and skip this if finished. Likely will want this to be more granular.

if [[ "${figuresFinished}" = "TRUE" ]]; then
  echo "Figure creation from prior run"
else
  if [ -d "${outdirectory}/Figures" ]; then
    rm -r ${outdirectory}/Figures
    mkdir ${outdirectory}/Figures
  else
    mkdir ${outdirectory}/Figures
  fi
  if [ -d "${outdirectory}/processed_tables" ]; then
    rm -r ${outdirectory}/processed_tables
    mkdir ${outdirectory}/processed_tables
  else
    mkdir ${outdirectory}/processed_tables
  fi
  
#Interpret sample metadata file for inputs:
highestgroupnum=0
rm -f ${workingdirectory}/${outdirectory}/chem_headers.txt

for ((f=1; f<=`awk '{print NF}' ${workingdirectory}/${outdirectory}/sample_metadata.txt | sort -nu | tail -n 1`; f++))
  do cat ${workingdirectory}/${outdirectory}/sample_metadata.txt | cut -f${f} > ${workingdirectory}/${outdirectory}/temp
    header=`head -n 1 ${workingdirectory}/${outdirectory}/temp`
    if [[ "$header" = "Sample" || "$header" = "lat" || "$header" = "long" ]]; then
      continue
      elif [[ "$header" = "replicates" ]]; then
      replicates="TRUE"
      elif [[ "$header" = "sites" ]]; then
      sites="TRUE"
      elif [[ "$header" = "controls" ]]; then
      controlspresent="TRUE"
      positiveCount=`cat ${workingdirectory}/${outdirectory}/temp | grep -c "positive"`
      negativeCount=`cat ${workingdirectory}/${outdirectory}/temp | grep -c "negative"`
      if [[ $positiveCount -ge 1 ]]; then
        controlPos="TRUE"
      fi
      if [[ $negativeCount -ge 1 ]]; then
        controlNeg="TRUE"
      fi
      elif [[ "$header" =~ "group" ]]; then
      groupsDefinedFlag="TRUE"
      numGrp=`echo $header | sed -E 's/group//'`
      if [[ $numGrp -gt $highestgroupnum ]]; then
        highestgroupnum=$numGrp
      fi
    else
      chemData="TRUE"
      echo $header >> ${workingdirectory}/${outdirectory}/chem_headers.txt
      locationChemHeaders="${workingdirectory}/${outdirectory}/chem_headers.txt"
    fi
  done
  rm ${workingdirectory}/${outdirectory}/temp
  numberGroupsDefined=$highestgroupnum

#echo "Replicates: $replicates"
#echo "Sites: $sites"
#echo "Controls: $controlspresent"
#echo "Positive Control: $controlPos"
#echo "Negative Control: $controlNeg"
#echo "Groups called: $groupsDefinedFlag"
#echo "Number of groups: $numberGroupsDefined"
#echo "Chem data: $chemData"
#echo "Location: $locationChemHeaders"

#Maps
mkdir ${outdirectory}/Figures/01_Maps
Rscript --vanilla ${metapipedir}/assets/maps.R ${workingdirectory}/${outdirectory}/Figures/01_Maps ${workingdirectory}/${outdirectory}/sample_metadata_forR.txt $replicates $sites ${workingdirectory}/${outdirectory}/ASV2Taxonomy/${outdirectory}_NO_UNKNOWNS_barchart.txt $filterPercent \
    1>> ${workingdirectory}/${outdirectory}/Figures/01_Maps/maps_rscript_out.log 2>&1

rm -f ${workingdirectory}/${outdirectory}/Figures/01_Maps/Rplot*

#Tables
perl ${metapipedir}/assets/barchart_filterLowAbund.pl -i ${outdirectory}/ASV2Taxonomy/${outdirectory}_barchart_forR.txt -f $filterPercent > ${outdirectory}/ASV2Taxonomy/${outdirectory}_barchart_forR_filtLowAbund_zzOther.txt

Rscript --vanilla ${metapipedir}/assets/process_tables.R ${workingdirectory}/${outdirectory} ${workingdirectory}/${outdirectory}/ASV2Taxonomy/ASVs_counts_NOUNKNOWNS.tsv ${workingdirectory}/${outdirectory}/ASV2Taxonomy/${outdirectory}_asvTaxonomyTable_NOUNKNOWNS.txt ${workingdirectory}/${outdirectory}/sample_metadata_forR.txt $filterPercent $controlspresent $filterLowQualSamples $filterPercentLowQualSamples ${workingdirectory}/${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy_NOUNKNOWNS.tsv $sites $replicates \
  1>> ${workingdirectory}/${outdirectory}/processed_tables/table_rscript_out.log 2>&1
perl ${metapipedir}/assets/filter_lowabundance_taxa.pl -a ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv -t ${outdirectory}/ASV2Taxonomy/${outdirectory}_asvTaxonomyTable_NOUNKNOWNS.txt -p $filterPercent > ${outdirectory}/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt

if [[ "${controlspresent}" = "TRUE" ]]; then
    if [[ "${filterLowQualSamples}" = "TRUE" ]]; then
      compareproc1=`cat ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_controlsRemoved.tsv`
      compareproc2=`cat ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv`
      compareproc3=`cat ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_controlsRemoved.tsv`
      compareproc4=`cat ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv`
    
      if [[ "${compareproc1}" = "${compareproc2}" ]]; then
        echo "Low quality sample filtering ON, but none met criteria to be removed"
        filterLowQualSamples=FALSE
        rm ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv
      fi
      if [[ "${compareproc3}" = "${compareproc4}" ]]; then
        rm ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv
      fi
    fi
else #CONTROLSPRESENT = FALSE
    if [[ "${filterLowQualSamples}" = "TRUE" ]]; then
      compareproc1=`cat ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy_NOUNKNOWNS.tsv`
      compareproc2=`cat ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv`
      compareproc3=`cat ${outdirectory}/ASV2Taxonomy/ASVs_counts_NOUNKNOWNS.tsv`
      compareproc4=`cat ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv`
    
      if [[ "${compareproc1}" = "${compareproc2}" ]]; then
        echo "Low quality sample filtering ON, but none met criteria to be removed"
        filterLowQualSamples=FALSE
        rm ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv
      fi
      if [[ "${compareproc3}" = "${compareproc4}" ]]; then
        rm ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv
      fi
    fi
fi

for f in ${workingdirectory}/${outdirectory}/processed_tables/sample_metadata*; do cat $f | sed -E 's/, /_/g' > ${f}_mod; mv ${f}_mod $f; done

#Phyloseq figures
mkdir -p ${outdirectory}/Figures/02_Barcharts/read_count
mkdir -p ${outdirectory}/Figures/02_Barcharts/relative_abundance
mkdir -p ${outdirectory}/Figures/03_Heatmaps/ASV_based
mkdir -p ${outdirectory}/Figures/03_Heatmaps/Taxonomy_merge_based
mkdir -p ${outdirectory}/Figures/04_Alpha_Diversity/ASV_based
mkdir -p ${outdirectory}/Figures/04_Alpha_Diversity/Taxonomy_merge_based
mkdir -p ${outdirectory}/Figures/05_Ordination/ASV_based/read_count
mkdir -p ${outdirectory}/Figures/05_Ordination/ASV_based/relative_abundance
mkdir -p ${outdirectory}/Figures/05_Ordination/Taxonomy_merge_based/read_count
mkdir -p ${outdirectory}/Figures/05_Ordination/Taxonomy_merge_based/relative_abundance
mkdir -p ${outdirectory}/Figures/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/read_count
mkdir -p ${outdirectory}/Figures/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/relative_abundance
mkdir -p ${outdirectory}/Figures/06_Network/ASV_based/read_count
mkdir -p ${outdirectory}/Figures/06_Network/ASV_based/relative_abundance
mkdir -p ${outdirectory}/Figures/06_Network/Taxonomy_merge_based/read_count
mkdir -p ${outdirectory}/Figures/06_Network/Taxonomy_merge_based/relative_abundance
mkdir -p ${outdirectory}/Figures/07_Rarefaction_Curves
mkdir -p ${outdirectory}/Figures/08_Taxa_of_interest/02_Barcharts/read_count
mkdir -p ${outdirectory}/Figures/08_Taxa_of_interest/02_Barcharts/relative_abundance
mkdir -p ${outdirectory}/Figures/08_Taxa_of_interest/03_Heatmaps/ASV_based
mkdir -p ${outdirectory}/Figures/08_Taxa_of_interest/03_Heatmaps/Taxonomy_merge_based
mkdir -p ${outdirectory}/Figures/08_Taxa_of_interest/06_Network/ASV_based/read_count
mkdir -p ${outdirectory}/Figures/08_Taxa_of_interest/06_Network/ASV_based/relative_abundance
mkdir -p ${outdirectory}/Figures/08_Taxa_of_interest/06_Network/Taxonomy_merge_based/read_count
mkdir -p ${outdirectory}/Figures/08_Taxa_of_interest/06_Network/Taxonomy_merge_based/relative_abundance

Rscript --vanilla ${metapipedir}/assets/phyloseq.R ${workingdirectory}/${outdirectory}/Figures ${workingdirectory} ${outdirectory} $controlspresent $filterLowQualSamples $replicates $sites $filterPercent $removeNA $providedTaxaOfInterest $groupsDefinedFlag $numberGroupsDefined $taxaOfInterestLevel $taxaOfInterestFile $chemData $locationChemHeaders \
  1>> ${workingdirectory}/${outdirectory}/Figures/phyloseq_rscript_out.log 2>&1

rm -f ${workingdirectory}/${outdirectory}/Figures/02_Barcharts/read_count/Rplots.pdf

Rscript --vanilla ${metapipedir}/assets/barchart_terminaltaxa.R ${workingdirectory}/${outdirectory}/Figures ${workingdirectory}/${outdirectory}/ASV2Taxonomy/${outdirectory}_barchart_forR.txt ${workingdirectory}/${outdirectory}/sample_order.txt ${workingdirectory}/${outdirectory}/ASV2Taxonomy/${outdirectory}_barchart_forR_filtLowAbund_zzOther.txt \
 1>> ${workingdirectory}/${outdirectory}/Figures/phyloseq_rscript_out.log 2>&1

rm -f ${workingdirectory}/${outdirectory}/Figures/02_Barcharts/relative_abundance/Rplots.pdf

#Replicate presence/absence abundance (presence in multiple replicates as sign of abundance)
#
#if [[ "${replicates}" = "TRUE" ]]; then
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#fi
#





  
  
fi #Final fi for if Figures folder present statement


echo "YOU MADE IT!"









