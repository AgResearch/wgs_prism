#!/bin/bash
#
# this prism supports a basic production and q/c wgs analysis, of data that is assumed to be 
# generated and hosted by AgResearch  - i.e. there are local dependencies 
# 
#

declare -a files_array

function get_opts() {

   DRY_RUN=no
   DEBUG=no
   HPC_TYPE=slurm
   FILES=""
   OUT_ROOT=""
   FORCE=no
   ANALYSIS=all

   help_text="
usage :\n 
./ag_wgs_qc_prism.sh  [-h] [-n] [-d] [-f] [-C hpctype] [-s SAMPLESHEET] [-a bcl2fastq|bcl2fastq_custom|fastqc|fasta_sample|fastq_sample|kmer_analysis|blast_analysis|all|html|common_sequence|clean] -O outdir  project1 project2 . . . \n
example:\n
./ag_wgs_qc_prism.sh -f -a bcl2fastq -O /dataset/gseq_processing/scratch/illumina/hiseq/190524_D00390_0462_ACDN6VANXX  -r 190524_D00390_0462_ACDN6VANXX -s /dataset/hiseq/active/190524_D00390_0462_ACDN6VANXX/SampleSheet.csv \n
./ag_wgs_qc_prism.sh -f -a fastqc -O /dataset/gseq_processing/scratch/illumina/hiseq/190524_D00390_0462_ACDN6VANXX  -r 190524_D00390_0462_ACDN6VANXX -s /dataset/hiseq/active/190524_D00390_0462_ACDN6VANXX/SampleSheet.csv  OLW_Plate1 \n
./ag_wgs_qc_prism.sh -a html  -O /dataset/gseq_processing/scratch/illumina/hiseq/190510_D00390_0457_AHYFCVBCX2 -r 190510_D00390_0457_AHYFCVBCX2 /dataset/hiseq/active/190510_D00390_0457_AHYFCVBCX2/SampleSheet.csv\n
"
   while getopts ":nhfO:C:r:a:s:" opt; do
   case $opt in
       n)
         DRY_RUN=yes
         ;;
       d)
         DEBUG=yes
         ;;
       f)
         FORCE=yes
         ;;
       h)
         echo -e $help_text
         exit 0
         ;;
       r)
         RUN=$OPTARG
         ;;
       s)
         SAMPLESHEET=$OPTARG
         ;;
       a)
         ANALYSIS=$OPTARG
         ;;
       C)
         HPC_TYPE=$OPTARG
         ;;
       O)
         OUT_ROOT=$OPTARG
         ;;
       \?)
         echo "Invalid option: -$OPTARG" >&2
         exit 1
         ;;
       :)
         echo "Option -$OPTARG requires an argument." >&2
         exit 1
         ;;
     esac
   done

   shift $((OPTIND-1))

   PROJECT_STRING=$@
   
   # when we run bcl2fastq there are not any projects yet
   if [ -z "$PROJECT_STRING" ]; then
      PROJECT_STRING="all"
   fi

   # this is needed because of the way we process args a "$@" - which
   # is needed in order to parse parameter sets to be passed to the
   # aligner (which are space-separated)
   declare -a projects="(${PROJECT_STRING})";
   NUM_PROJECTS=${#projects[*]}
   for ((i=0;$i<$NUM_PROJECTS;i=$i+1)) do
      projects_array[$i]=${projects[$i]}
   done

}


function check_opts() {
   if [ -z "$SEQ_PRISMS_BIN" ]; then
      echo "please set SEQ_PRISMS_BIN environment variable"
      exit 1
   fi

   if [ -z "$WGS_PRISM_BIN" ]; then
      echo "please set WGS_PRISM_BIN environment variable"
      exit 1
   fi

   if [ ! -d $OUT_ROOT ]; then
      echo "out_dir $OUT_ROOT not found"
      exit 1
   fi

   if [ ! -f $SAMPLESHEET ]; then
      echo "$SAMPLESHEET not found"
      exit 1
   fi


   if [[ $HPC_TYPE != "local" && $HPC_TYPE != "slurm" ]]; then
      echo "HPC_TYPE must be one of local, slurm"
      exit 1
   fi

   if [[ ( $ANALYSIS != "all" ) && ( $ANALYSIS != "bcl2fastq" ) && ( $ANALYSIS != "fastqc" ) && ( $ANALYSIS != "bcl2fastq_custom" ) && ( $ANALYSIS != "clean" ) && ( $ANALYSIS != "fasta_sample" ) && ( $ANALYSIS != "kmer_analysis" ) && ( $ANALYSIS != "blast_analysis" ) && ( $ANALYSIS != "annotation" )  && ( $ANALYSIS != "html" ) && ( $ANALYSIS != "clientreport" )  && ( $ANALYSIS != "fastq_sample" ) && ( $ANALYSIS != "common_sequence" ) ]]; then
      echo "analysis must be one of all, demultiplex, kgd , unblind, kmer_analysis, allkmer_analysis, blast_analysis , common_sequencs, clean "
      exit 1
   fi

}

function echo_opts() {
  echo OUT_ROOT=$OUT_ROOT
  echo DRY_RUN=$DRY_RUN
  echo DEBUG=$DEBUG
  echo HPC_TYPE=$HPC_TYPE
  echo ANALYSIS=$ANALYSIS

}


#
# edit this method to set required environment (or set up
# before running this script)
#
function configure_env() {
   export CONDA_ENVS_PATH=$CONDA_ENVS_PATH:/dataset/bioinformatics_dev/active/conda-env
   cd $WGS_PRISM_BIN
   cp ag_wgs_qc_prism.sh $OUT_ROOT
   cp ag_wgs_qc_prism.mk $OUT_ROOT

   echo "
max_tasks=50
jobtemplatefile = \"$WGS_PRISM_BIN/etc/wgs_qc_slurm_array_job\"
" > $OUT_ROOT/tardis.toml

   echo "
conda activate bifo-essential
" > $OUT_ROOT/bifo-essential_env.inc

   echo "
export CONDA_ENVS_PATH=$CONDA_ENVS_PATH
conda activate bioconductor
" > $OUT_ROOT/configure_bioconductor_env.src


   cd $OUT_ROOT
}


function check_env() {
   if [ -z "$SEQ_PRISMS_BIN" ]; then
      echo "SEQ_PRISMS_BIN not set - exiting"
      exit 1
   fi
   if [ -z "$WGS_PRISM_BIN" ]; then
      echo "WGS_PRISM_BIN not set - exiting"
      exit 1
   fi

}


function get_targets() {
   # make target monikers  and write associated
   # wrapper, which will be called by make

   rm -f $OUT_ROOT/*_targets.txt

   for ((j=0;$j<$NUM_PROJECTS;j=$j+1)) do
      projectname=${projects_array[$j]}
      samplesheet_base=`basename $SAMPLESHEET .csv`

      if [ $ANALYSIS != "bcl2fastq" ]; then
         project_moniker=${RUN}.${projectname}.$samplesheet_base
      else
         project_moniker=${RUN}.$samplesheet_base
      fi

      for analysis_type in all bcl2fastq fastqc clean kmer_analysis blast_analysis fasta_sample fastq_sample annotation common_sequence; do
         echo $OUT_ROOT/$project_moniker.$analysis_type  >> $OUT_ROOT/${analysis_type}_targets.txt
         script=$OUT_ROOT/${project_moniker}.${analysis_type}.sh
         if [ -f $script ]; then
            if [ ! $FORCE == yes ]; then
               echo "found existing wg script $script  - will re-use (use -f to force rebuild of scripts) "
               continue
            fi
         fi
      done

      ############### bcl2fastq script 
      echo "#!/bin/bash
cd $OUT_ROOT
mkdir -p $samplesheet_base
$SEQ_PRISMS_BIN/sequencing_qc_prism.sh -a bcl2fastq -O $OUT_ROOT/$samplesheet_base $SAMPLESHEET > $OUT_ROOT/$samplesheet_base/bcl2fastq.log  2>&1
if [ \$? != 0 ]; then
   echo \"warning bcl2fastq of $samplesheet_base returned an error code\"
      exit 1
   fi
      " > $OUT_ROOT/${project_moniker}.bcl2fastq.sh
      chmod +x $OUT_ROOT/${project_moniker}.bcl2fastq.sh


      ############### bcl2fastq_custom script
      echo "#!/bin/bash
cd $OUT_ROOT
mkdir -p $samplesheet_base
$SEQ_PRISMS_BIN/sequencing_qc_prism.sh -a bcl2fastq -B \"--ignore-missing-bcls\"  -O $OUT_ROOT/$samplesheet_base $SAMPLESHEET > $OUT_ROOT/$samplesheet_base/bcl2fastq.log  2>&1
if [ \$? != 0 ]; then
   echo \"warning custom bcl2fastq of $samplesheet_base returned an error code\"
      exit 1
   fi
      " > $OUT_ROOT/${project_moniker}.bcl2fastq_custom.sh
      chmod +x $OUT_ROOT/${project_moniker}.bcl2fastq_custom.sh


      ############### fastqc script 
      echo "#!/bin/bash
cd $OUT_ROOT
mkdir -p $samplesheet_base/$projectname

# collate a combined R1 and R2 file and run fastqc on those
for pair_moniker in R2 R1 ; do
   #ls $OUT_ROOT/$samplesheet_base/bcl2fastq/$projectname/*_\${pair_moniker}_*.fastq.gz | grep -v Blank_ > $OUT_ROOT/$samplesheet_base/$projectname/\${pair_moniker}.list
   find $OUT_ROOT/$samplesheet_base/bcl2fastq/$projectname -name \"*_\${pair_moniker}_*.fastq.gz\" | grep -v Blank_ > $OUT_ROOT/$samplesheet_base/$projectname/\${pair_moniker}.list
   rm $OUT_ROOT/$samplesheet_base/$projectname/all_demultiplexed\${pair_moniker}.fastq
   for file in \`cat $OUT_ROOT/$samplesheet_base/$projectname/\${pair_moniker}.list\`; do
      gunzip -c \$file >> $OUT_ROOT/$samplesheet_base/$projectname/all_demultiplexed\${pair_moniker}.fastq
   done
   $SEQ_PRISMS_BIN/sequencing_qc_prism.sh -C local -a fastqc -O $OUT_ROOT/$samplesheet_base/$projectname $OUT_ROOT/$samplesheet_base/$projectname/all_demultiplexed\${pair_moniker}.fastq >> $OUT_ROOT/$samplesheet_base/$projectname/fastqc.log 2>&1
   #rm $OUT_ROOT/$samplesheet_base/all_demultiplexed${pair_moniker}.fastq
done


# also run individual fastqc on each sample
$SEQ_PRISMS_BIN/sequencing_qc_prism.sh -a fastqc -O $OUT_ROOT/$samplesheet_base/$projectname \`cat $OUT_ROOT/$samplesheet_base/$projectname/*.list\`  >  $OUT_ROOT/$samplesheet_base/$projectname/fastqc_by_sample.log

if [ \$? != 0 ]; then
   echo \"warning fastqc of $projectname returned an error code\"
   exit 1
fi
     " > $OUT_ROOT/${project_moniker}.fastqc.sh
     chmod +x $OUT_ROOT/${project_moniker}.fastqc.sh


     ################ fastq_sample script 
     echo "#!/bin/bash
cd $OUT_ROOT
mkdir -p $samplesheet_base/$projectname/sample

find $OUT_ROOT/$samplesheet_base/bcl2fastq/$projectname -name \"*.fastq.gz\"  | grep -v \".005\" > $OUT_ROOT/$samplesheet_base/$projectname/fastq_sample.list

$SEQ_PRISMS_BIN/sample_prism.sh -s .005 -M 50000 -a fastq -O $OUT_ROOT/$samplesheet_base/$projectname/sample \`cat $OUT_ROOT/$samplesheet_base/$projectname/fastq_sample.list \`


if [ \$? != 0 ]; then
   echo \"warning sampling of seqs in $OUT_ROOT/$samplesheet_base/$projectname/fastq_sample.list returned an error code\"
   exit 1
fi
     " >  $OUT_ROOT/${project_moniker}.fastq_sample.sh
      chmod +x $OUT_ROOT/${project_moniker}.fastq_sample.sh

     ################ fasta_sample script
     echo "#!/bin/bash
cd $OUT_ROOT
mkdir -p $samplesheet_base/$projectname/sample

find $OUT_ROOT/$samplesheet_base/bcl2fastq/$projectname -name \"*.fastq.gz\"  | grep -v \".005\" > $OUT_ROOT/$samplesheet_base/$projectname/fasta_sample.list

$SEQ_PRISMS_BIN/sample_prism.sh -s .0001 -M 1000 -a fasta -O $OUT_ROOT/$samplesheet_base/$projectname/sample \`cat $OUT_ROOT/$samplesheet_base/$projectname/fasta_sample.list \`

if [ \$? != 0 ]; then
   echo \"warning sampling of seqs in $OUT_ROOT/$samplesheet_base/$projectname/fasta_sample.list returned an error code\"
   exit 1
fi
     " >  $OUT_ROOT/${project_moniker}.fasta_sample.sh
      chmod +x $OUT_ROOT/${project_moniker}.fasta_sample.sh

     ################ kmer summary script
     echo "#!/bin/bash
cd $OUT_ROOT
mkdir -p $samplesheet_base/$projectname/kmer_analysis
$SEQ_PRISMS_BIN/kmer_prism.sh -C $HPC_TYPE  -a fasta -p \"-k 6 -A\" -O $OUT_ROOT/$samplesheet_base/$projectname/kmer_analysis $OUT_ROOT/$samplesheet_base/$projectname/sample/*.fasta.gz 
if [ \$? != 0 ]; then
   echo \"warning, kmer analysis of $OUT_ROOT/$samplesheet_base/$projectname/sample/*.fastq.gz  returned an error code\"
   exit 1
fi
     " >  $OUT_ROOT/${project_moniker}.kmer_analysis.sh
      chmod +x $OUT_ROOT/${project_moniker}.kmer_analysis.sh


     ################ clean script
     echo "#!/bin/bash
cd $OUT_ROOT
rm -rf $OUT_ROOT/$samplesheet_base/$projectname
rm -f *.${samplesheet_base}.*
     " >  $OUT_ROOT/${project_moniker}.clean.sh
      chmod +x $OUT_ROOT/${project_moniker}.clean.sh



     ################ blast script 
### note ###
# Re "Misunderstood parameter of NCBI BLAST impacts the correctness of bioinformatics workflows", Nidhi Shah  Michael G Nute  Tandy Warnow  Mihai Pop
# and our use of "-max_target_seqs 1"
# In the present context, the top hit returned from each (randomly sampled) sequence, from each sequenced biological sample, 
# is used to prepare a numeric profile vector for each file, with the semantic details of the hits discarded.
# The numeric vectors are then input to unsupervised machine learning - for example clustered
# - so that we can highlight how similar or dissimilar new files are to previous files, and to each other.
# It is not necessary for our purpose here that the hit returned , is the best (i.e. lowest evalue) in the database.
# (This ("non-semantic") approach does depend on the same database version being used throughout
# the series of files - and this would be true even if this blast parameter behaved as intuitively 
# expected - i.e. returned the actual best hit in the database).
############
     echo "#!/bin/bash
cd $OUT_ROOT
mkdir -p $samplesheet_base/$projectname/blast
# configure a custom slurm batch job that will specify medium memory 
cp $WGS_PRISM_BIN/etc/medium_mem_slurm_array_job $OUT_ROOT
echo \"
jobtemplatefile = \\\"$OUT_ROOT/medium_mem_slurm_array_job\\\"
max_tasks = 60
\" > $OUT_ROOT/$samplesheet_base/$projectname/blast/tardis.toml
# run blast
$SEQ_PRISMS_BIN/align_prism.sh -C $HPC_TYPE -m 60 -a blastn -r nt -p \"-evalue 1.0e-10  -dust \\'20 64 1\\' -max_target_seqs 1 -outfmt \\'7 qseqid sseqid pident evalue staxids sscinames scomnames sskingdoms stitle\\'\" -O $OUT_ROOT/$samplesheet_base/$projectname/blast $OUT_ROOT/$samplesheet_base/$projectname/fasta_small_lowdepthsample/*.fasta
if [ \$? != 0 ]; then
   echo \"warning , blast  of $OUT_ROOT/$samplesheet_base/$projectname/fasta_small_lowdepthsample returned an error code\"
   exit 1
fi
     " >  $OUT_ROOT/${project_moniker}.blast_analysis.sh 
      chmod +x $OUT_ROOT/${project_moniker}.blast_analysis.sh


     ################ annotation script 
     echo "#!/bin/bash
cd $OUT_ROOT
# summarise species from blast results 
$SEQ_PRISMS_BIN/annotation_prism.sh -C $HPC_TYPE -w tag_count -a taxonomy -O $OUT_ROOT/$samplesheet_base/$projectname/blast $OUT_ROOT/$samplesheet_base/$projectname/blast/*.results.gz  
return_code1=\$?
# summarise descriptions from blast results 
rm -f $OUT_ROOT/$samplesheet_base/$projectname/blast/*.annotation_prism
$SEQ_PRISMS_BIN/annotation_prism.sh -C $HPC_TYPE -w tag_count -a description -O $OUT_ROOT/$samplesheet_base/$projectname/blast $OUT_ROOT/$samplesheet_base/$projectname/blast/*.results.gz  
return_code2=\$?
# provide unblinded frequency tables
for file in  $OUT_ROOT/$samplesheet_base/$projectname/blast/frequency_table.txt $OUT_ROOT/$samplesheet_base/$projectname/blast/locus_freq.txt ; do
   if [ -f \$file ]; then
      if [ ! -f \$file.blinded ]; then
         cp -p \$file \$file.blinded
      fi
      cat \$file.blinded | sed -f $OUT_ROOT/${project_moniker}.unblind.sed > \$file
   fi
done
if [[ ( \$return_code1 != 0 ) || ( \$return_code2 != 0 ) ]]; then
   echo \"warning, summary of $OUT_ROOT/$samplesheet_base/$projectname/blast returned an error code\"
   exit 1
fi
     " >  $OUT_ROOT/${project_moniker}.annotation.sh 
      chmod +x $OUT_ROOT/${project_moniker}.annotation.sh 

   done
}



function fake_prism() {
   echo "dry run ! 

   "
   exit 0
}

function run_prism() {
   cd $OUT_ROOT

   make -f ag_wgs_qc_prism.mk -d -k  --no-builtin-rules -j 16 `cat $OUT_ROOT/${ANALYSIS}_targets.txt` > $OUT_ROOT/${ANALYSIS}.log 2>&1

   # run summaries
}

function html_prism() {
   mkdir -p $OUT_ROOT/html

   samplesheet_base=`basename $SAMPLESHEET .csv`

   for folder_path in $OUT_ROOT/$samplesheet_base/bcl2fastq/* ;  do
      if [ ! -d $folder_path ]; then
         continue
      fi
      projectname=`basename $folder_path`
      if [[ ( $projectname == "Reports" ) || ( $projectname == "Stats" ) ]]; then 
         continue
      fi 

      mkdir -p $OUT_ROOT/html/$projectname
      cp -pR $OUT_ROOT/$samplesheet_base/$projectname/fastqc $OUT_ROOT/html/$projectname
      
      mkdir -p  $OUT_ROOT/html/$projectname/kmer_analysis
      rm $OUT_ROOT/html/$projectname/kmer_analysis/*
      for file in $OUT_ROOT/$samplesheet_base/$projectname/kmer_analysis/*.jpg $OUT_ROOT/$samplesheet_base/$projectname/kmer_analysis/*.txt ; do
          cp -s $file $OUT_ROOT/html/$projectname/kmer_analysis
      done

      mkdir -p  $OUT_ROOT/html/$projectname/blast
      rm $OUT_ROOT/html/$projectname/blast/*
      for file in $OUT_ROOT/$samplesheet_base/$projectname/blast/*.jpg $OUT_ROOT/$samplesheet_base/$projectname/blast/taxonomy*clusters.txt $OUT_ROOT/$samplesheet_base/$projectname//blast/frequency_table.txt  $OUT_ROOT/$cohort/blast/locus_freq.txt; do
         cp -s $file $OUT_ROOT/html/$projectname/blast
      done

      mkdir -p $OUT_ROOT/html/$projectname/common_sequence
      rm $OUT_ROOT/html/$projectname/common_sequence/all_common_sequence.txt
      rm $OUT_ROOT/html/$projectname/common_sequence/preview_common_sequence.txt
      for file in $OUT_ROOT/$samplesheet_base/$projectname/kmer_analysis/*.k6A.log ; do
         base=`basename $file .k6A.log`
         echo "

         sample: $base" >> $OUT_ROOT/html/$projectname/common_sequence/preview_common_sequence.txt
         echo "

         sample: $base" >> $OUT_ROOT/html/$projectname/common_sequence/all_common_sequence.txt
         grep assembled_by_distinct $file | sed 's/assembled_by_distinct//g' >> $OUT_ROOT/html/$projectname/common_sequence/all_common_sequence.txt
         grep assembled_by_distinct $file | head -10 | sed 's/assembled_by_distinct//g' >> $OUT_ROOT/html/$projectname/common_sequence/preview_common_sequence.txt
      done
   done

   mkdir -p $OUT_ROOT/html/bcl2fastq
   cp -pR $OUT_ROOT/$samplesheet_base/bcl2fastq/Reports/html/* $OUT_ROOT/html/bcl2fastq

   # make peacock page which mashes up plots, output files etc.
   $WGS_PRISM_BIN/make_sample_pages.py -r $RUN -o $OUT_ROOT/html/peacock.html
}

function clientreport_prism() {
   if [ ! -d $OUT_ROOT/html ]; then 
      echo "could not find $OUT_ROOT/html (please generate the html summaries first)"
      exit 1
   fi

   $WGS_PRISM_BIN/make_clientsamplesheet_pages.py -r $RUN -o report.html
   for ((j=0;$j<$NUM_PROJECTS;j=$j+1)) do
      set -x
      SAMPLESHEET=${SAMPLESHEETs_array[$j]}
      cd  $OUT_ROOT/html/$SAMPLESHEET 
      rm -f report.zip report.tar*
      tar -cv --auto-compress --dereference -f report.tar.gz  --files-from=report.html.manifest
      cat report.html.manifest | zip -@ report
      set +x
   done
}

function clean() {
   echo "cleaning up tardis working folders..."
   nohup find $OUT_ROOT -name "tardis_*" -type d -exec rm -r {} \; > $OUT_ROOT/clean.log 2>&1 &
}


function main() {
   get_opts "$@"
   check_opts
   echo_opts
   check_env
   configure_env

   if [ $ANALYSIS == "html" ]; then
      html_prism
   elif [ $ANALYSIS == "clientreport" ]; then
      clientreport_prism
   else
      get_targets
      if [ $DRY_RUN != "no" ]; then
         fake_prism
      else
         run_prism
         if [ $? == 0 ] ; then
            clean
            echo "* done clean *"  # mainly to yield zero exit code
         else
            echo "error state from  run - skipping html page generation and clean-up"
            exit 1
         fi
      fi
   fi
}


set -x
main "$@"
set +x
