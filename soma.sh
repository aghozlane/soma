#!/bin/bash
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html
# ------------------------------------------------------------------
# Authors: Amine Ghozlane (amine.ghozlane@jouy.inra.fr)
#          Mathieu Almeida
#          Florian Plaza Oñate
# Title:  Assemble_SOLiD_lite-metagenome
# Description : Assembly pipeline for SOLID reads with annotation
# ------------------------------------------------------------------

function say_parameters {
    # echo blue
    echo -e "\e[34m# $1\e[0m" >&2
}

function say {
    # echo green
    echo -e "\033[1;32m* $1\033[0m" >&2
}

function error {
    # echo red
    echo -e  "\e[31m* $1\e[0m" >&2
}

function check_log {
    # check if log file is not empty
    if [ -s $1 ]
    then
        error "$1 is not empty !"
        exit 1
    fi
}

function check_integer {
    if [[ $1 != [0-9]* ]]
    then
        error "\"$1\" is not an integer value"
        exit 1
    fi
}

function check_file {
    # Check if result is well produced
    if [ ! -f $1 ] && [ ! -s $1 ]
    then
	error "File \"$1\" does not exist or is empty !"
	exit 1
    fi
}

function check_dir {
    # Check if directory doesnt exist
    if [ ! -d $1 ]
    then
        mkdir $1
        if [ ! -d $1 ]
        then
            error "The program cannot create the directory \"$1\""
            exit 1
        fi
    fi
}

display_help() {
    if [ "$1" -eq "0" ]
    then
        echo """Usage :
case 1 - Assembly: $0 -i <read.csfasta> -o </path/to/result/directory/> --assembly (--hsize 23 -r 20000000 for 35nt reads or --hsize 25 -r 100000000 for 50nt reads --abundance_filtering --lowerCount 3 --occurence 3 --kmer 15  for low quality sequencing)
case 2 - Taxonomic annotation/MOTU: $0 -g <gene.fasta> -o </path/to/result/directory/> --tax_annotation
case 3 - Functional annotation: $0 -p <protein.fasta> -o </path/to/result/directory/> --func_annotation
all : $0 -i <read.csfasta> -q <rea.quality> -o </path/to/result/directory/>
        """
        echo "Use '--help' to print detailed descriptions of command line arguments"
    else
        display_parameters
    fi
    exit
}

display_parameters() {
   # Display the parameters of the analysis
   say_parameters "Sample csfasta [-i] :"
   echo $input >&2
   say_parameters "Input gene [-g] :"
   echo $input_gene >&2
   say_parameters "Input protein [-p] :"
   echo $input_protein >&2
   say_parameters "Output directory [-o] :"
   echo $resultDir  >&2
   say_parameters "Filtering reads from the databases :"
   for db in ${filterRef[@]}
   do
    echo $db  >&2
   done
   say_parameters "Assembly parameters :"
   echo "Contig length threshold [-c] = $contigLengthThreshold" >&2
   echo "Read trim length [-t] = $readTrimLength" >&2
   echo "Number of mismatch for the mapping [-m] = $NbMismatchMapping">&2
   echo "Denovo2 options = $denovo2Options">&2
   echo "Reference length [-r] = $refLength">&2
   echo "Gene length threshold [-c] = $geneLengthThreshold" >&2
   echo "Number of process [-n|--NbProc] = $NbProc" >&2
   echo "Abundance filtering [--abundance_filtering] = $abundance_filtering"  >&2
   echo "Kmer size [--kmer] = $kmer"  >&2
   echo "Lower count [--lowerCount] = $lowerCount"  >&2
   echo "Occurence of the kmer in the read [--occurence]  = $occurence"  >&2
   say_parameters "Gene prediction :"
   echo "Use metagenemark [--metagenemark] = $metagenemark"
   echo "Gene type [-f] = $filterGene"
   say_parameters "-Functional annotation-"
   say_parameters "Egg/Nog annotation :"
   echo "Egg/Nog sequence database [--eggnogRefBlast] = $eggnogRefBlast" >&2
   echo "Egg/Nog description [--eggnogdescRef] = $eggnogdescRef" >&2
   echo "Egg/Nog members [--eggnogmemRef] = $eggnogmemRef" >&2
   echo "Egg/Nog functional category description [--eggnogfunccats] = $eggnogfunccats" >&2
   echo "Egg/Nog functional category [--eggnogfunccatRef] = $eggnogfunccatRef" >&2
   say_parameters "Cazy annotation :"
   echo "Cazy HMM db [--cazyRefHmmer] = $cazyRefHmmer" >&2
   echo "Cazy family information [--famInfoRef] = $famInfoRef" >&2
}

function timer()
{
   if [[ $# -eq 0 ]]; then
         echo $(date '+%s')
   else
      local  stime=$1
      etime=$(date '+%s')
      if [[ -z "$stime" ]]; then stime=$etime; fi
      dt=$((etime - stime))
      ds=$((dt % 60))
      dm=$(((dt / 60) % 60))
      dh=$((dt / 3600))
      printf '%d:%02d:%02d' $dh $dm $ds
  fi
}

# Assembly metagenome path
SCRIPTPATH=$(readlink -f $(dirname "${BASH_SOURCE[0]}"))
DATAPATH=$SCRIPTPATH/data/

#############
# Databases #
#############
filterRef=("$DATAPATH/Homo_sapiens_2014_02_04/homo_sapiens.fna" "$DATAPATH/barcodes/barcodes.fna" )
eggnogRefBlast="$DATAPATH/eggnog_v3/eggnog_v3.faa"
eggnogdescRef="$DATAPATH/eggnog_v3/info/description/"
eggnogmemRef="$DATAPATH/eggnog_v3/info/members/"
eggnogfunccatRef="$DATAPATH/eggnog_v3/info/funccat/"
eggnogfunccats="$DATAPATH/eggnog_v3/info/eggnogv3.funccats.txt"
cazyRefHmmer="$DATAPATH/cazy_2013_05_11/dbCAN-fam-HMMs.txt"
famInfoRef="$DATAPATH/cazy_2013_05_11/info/FamInfo.txt"
gitaxidnucl="$DATAPATH/ncbi_nt/info/gi_taxid_nucl.dmp"
nodes="$DATAPATH/ncbi_nt/info/nodes.dmp"
names="$DATAPATH/ncbi_nt/info/names.dmp"
ntBlast="$DATAPATH/ncbi_nt/nt"



#######################
# Assembly Parameters #
#######################
refLength=10000000
# Trim des 0 dernieres bases pour reads de 50 nt
readTrimLength=0
NbMismatchMapping=1
contigLengthThreshold=300
NbProc=$(grep -c ^processor /proc/cpuinfo)
geneLengthThreshold=60
# Gene filtering
filterGene=11
# Activity
all=1
assembly=0
tax_annotation=0
func_annotation=0
maxTargetSeqs=20
evalueTaxAnnot="1e-3"
evalueFuncAnnot="1e-5"
numberBestannotation=1
denovo2Options=""
abundance_filtering=0
lowerCount=3
kmer=10
occurence=1
metagenemark=0


############
# Programs #
############
FilterFasta="$SCRIPTPATH/FilterFasta/FilterFasta.py"
extractimomi="$SCRIPTPATH/ExtractIMOMIAnnotation/ExtractIMOMIAnnotation.py"
trimReads="$SCRIPTPATH/TrimReads/trimReads.py"
filterProdigal="$SCRIPTPATH/FilterProdigal/FilterProdigal.py"
extractEggNog="$SCRIPTPATH/ExtractEggNog/ExtractEggNog.py"
extractCazy="$SCRIPTPATH/ExtractCazy/ExtractCazy.py"
extractProteins="$SCRIPTPATH/ExtractProteins/ExtractProteins.py"
extractKegg="$SCRIPTPATH/ExtractKegg/ExtractKegg.py"
extractNCBIDB="$SCRIPTPATH/ExtractNCBIDB/ExtractNCBIDB.py"
grabcataloguesequence="$SCRIPTPATH/grab_catalogue_sequence/grab_catalogue_sequence.py"
#gettaxonomy="$SCRIPTPATH/get_taxonomy/get_taxonomy.py"
gettaxonomy="$SCRIPTPATH/get_taxonomy/get_taxonomy"
extractMetagenemark="$SCRIPTPATH/ExtractMetagenemark/ExtractMetagenemark.py"
# git clone https://github.com/fplaza/filter-reads.git
filterReads="$SCRIPTPATH/filter-reads/filter-reads"
# git clone https://github.com/fplaza/csfasta-fasta-conv.git
csfastatofasta="$SCRIPTPATH/csfasta-fasta-conv/csfasta-fasta-conv"
# Assembly
export denovo2="$SCRIPTPATH/denovo2/"
# Blast
blastp="$SCRIPTPATH/ncbi-blast-2.2.31/blastp" #$(which blastp) 
blastn="$SCRIPTPATH/ncbi-blast-2.2.31/blastn" #$(which blastn)
# Bowtie
bowtie="$SCRIPTPATH/bowtie-1.1.2/bowtie" #$(which bowtie)
bowtie_build="$SCRIPTPATH/bowtie-1.1.2/bowtie-build" #$(which bowtie-build) 
# Prodigal
prodigal="$SCRIPTPATH/Prodigal-2.6.2/prodigal" # $(which prodigal)
# Hmmer
hmmscan="$SCRIPTPATH/hmmer-3.1b2/hmmscan" #$(which hmmscan)
# Jellyfish
jellyfish="$SCRIPTPATH/jellyfish/bin/jellyfish"
# Metagenemark
gmhmmp="$SCRIPTPATH/MetaGeneMark/mgm/gmhmmp"
metagenemarkmod="$SCRIPTPATH/MetaGeneMark/mgm/MetaGeneMark_v1.mod"
# FetchMG
fetchmg="$SCRIPTPATH/fetchMG/fetchMG.pl"

########
# Main #
########

# Execute getopt on the arguments passed to this program, identified by the special character $@
PARSED_OPTIONS=$(getopt -n "$0"  -o hi:o:n:m:t:c:r:l:f:p:g:q: --long "help,input:,quality:,input_protein:,input_gene:,output:,NbProc:,NbMismatchMapping:,readTrimLength:,contigLengthThreshold:,refLength:,geneLengthThreshold:,filterGene:,metagenemark,hsize:,NoCorrection,evalueTaxAnnot:,evalueFuncAnnot:,numberBestannotation:,maxTargetSeqs:,referenceGenome:,assembly,tax_annotation,func_annotation,abundance_filtering,kmer:,lowerCount:,occurence:,imomiRefBlast:,imomiIDtoGenome:,markerCOGRef:,markerCOGRefBlast:,eggnogRef:,eggnogRefBlast:,eggnogdescRef:,eggnogmemRef:,eggnogfunccats:,eggnogfunccatRef:,cazyRefHmmer:,famInfoRef:,keggRef:,keggRefBlast:,koRef:,koGeneKoRef:,koGeneUniprotRef:,koUniprotRef:,koGeneGeneidRef:,koGeneGiRef:,koGenePathwayRef:"  -- "$@")

#Check arguments
if [ $# -eq 0 ]
then
    display_help 0
fi

#Bad arguments, something has gone wrong with the getopt command.
if [ $? -ne 0 ];
then
    display_help 1
    exit 1
fi

# A little magic, necessary when using getopt.
eval set -- "$PARSED_OPTIONS"

# Get Cmd line arguments depending on options
while true;
do
  case "$1" in
    -h)
        display_help 0
        shift;;
    --help)
        display_help 1
        shift;;
    -i|--input)
        check_file $2
        input=$2
        shift 2;;
    -q|--quality)
        check_file $2
        quality=$2
        shift 2;;
    -p|--input_protein)
        check_file $2
        input_protein=$2
        shift 2;;
    -g|--input_gene)
        check_file $2
        input_gene=$2
        shift 2;;
    -o|--output)
        resultDir=$2
        logDir=$resultDir/log/
        errorlogDir=$resultDir/error_log/
        check_dir $resultDir
        check_dir $logDir
        check_dir $errorlogDir
        shift 2;;
    -n|--NbProc)
        check_integer $2
        NbProc=$2
        shift 2;;
    -m|--NbMismatchMapping)
        check_integer $2
        NbMismatchMapping=$2
        shift 2;;
    -t|--readTrimLength)
        check_integer $2
        readTrimLength=$2
        shift 2;;
    -c|--contigLengthThreshold)
        check_integer $2
        contigLengthThreshold=$2
        shift 2;;
    -r|--refLength)
        check_integer $2
        refLength=$2
        shift 2;;
    -l|--geneLengthThreshold)
        check_integer $2
        geneLengthThreshold=$2
        shift 2;;
    -f|--filterGene)
        check_integer $2
        filterGene=$2
        shift 2;;
    --metagenemark)
        metagenemark=1
        shift;;
    --hsize)
        check_integer $2
        denovo2Options+="-hsize $2 "
        shift 2;;
    --NoCorrection)
        denovo2Options+="-NO_CORRECTION "
        shift;;
    --evalueTaxAnnot)
        evalueTaxAnnot=$2
        shift 2;;
    --evalueFuncAnnot)
        evalueFuncAnnot=$2
        shift 2;;
    --numberBestannotation)
        check_integer $2
        numberBestannotation=$2
        shift 2;;
    --maxTargetSeqs)
        check_integer $2
        maxTargetSeqs=$2
        shift 2;;
    --referenceGenome)
        check_file $2
        referenceGenome=$2
        shift 2;;
    --assembly)
        assembly=1
        shift;;
    --tax_annotation)
        tax_annotation=1
        shift;;
    --func_annotation)
        func_annotation=1
        shift;;
    --abundance_filtering)
        abundance_filtering=1
        shift;;
    --kmer)
        check_integer $2
        kmer=$2
        shift 2;;
    --lowerCount)
        check_integer $2
        lowerCount=$2
        shift 2;;
    --occurence)
        check_integer $2
        occurence=$2
        shift 2;;
    --imomiRefBlast)
        imomiRefBlast=$2
        shift 2;;
    --imomiIDtoGenome)
        check_file $2
        imomiIDtoGenome=$2
        shift 2;;
    --markerCOGRef)
        check_file $2
        markerCOGRef=$2
        shift 2;;
    --markerCOGRefBlast)
        markerCOGRefBlast=$2
        shift 2;;
    --eggnogRef)
        check_file $2
        eggnogRef=$2
        shift 2;;
    --eggnogRefBlast)
       eggnogRefBlast=$2
       shift 2;;
    --eggnogdescRef)
        check_file $2
        eggnogdescRef=$2
        shift 2;;
    --eggnogmemRef)
        check_file $2
        eggnogmemRef=$2
        shift 2;;
    --eggnogfunccats)
        check_file $2
        eggnogfunccats=$2
        shift 2;;
    --eggnogfunccatRef)
        check_file $2
        eggnogfunccatRef=$2
        shift 2;;
    --cazyRefHmmer)
        check_file $2
        cazyRefHmmer=$2
        shift 2;;
    --famInfoRef)
        check_file $2
        famInfoRef=$2
        shift 2;;
    --)
      shift
      break;;
  esac
done

# Detect activity
if [ "$assembly" -eq "1" ] || [ "$tax_annotation" -eq "1" ] || [ "$func_annotation" -eq "1" ]
then
    all=0
fi

if [ "$resultDir" = "" ]
then
    error "Please indicate the output directory."
    exit 1
fi

# Check sample
if [ -f "$input" ]
then
    filename=$(basename "$input")
    extension=".${filename##*.}"
    if [ "$extension" != ".csfasta" ]
    then
        error "The input file should be a csfasta file."
        display_help
    fi
    SamplePath=$(dirname $input)
    SampleName="${filename%.*}"
    #SampleName=${filename::-8}
    # Set of gene and protein data
    input_gene=${resultDir}/${SampleName}_gene_${geneLengthThreshold}_filtered.fna
    input_protein=${resultDir}/${SampleName}_prot_filtered.faa
elif [ -f "$input_gene" ]
then
    filename=$(basename "$input_gene")
    SamplePath=$(dirname $input_gene)
    SampleName="${filename%.*}"
    if [ ! -f "$input_protein" ]
    then
        input_protein=${resultDir}/${SampleName}_prot_filtered.faa
    fi
elif [ -f "$input_protein" ]
then
    filename=$(basename "$input_protein")
    SamplePath=$(dirname $input_protein)
    SampleName="${filename%.*}"
else
    error "No input file !"
    exit 1
fi

# display parameters
display_parameters

# Start timer
say "Start analysis"
wall_time=$(timer)

# Assembly
if [ "$assembly" -eq "1" ] || [ "$all" -eq "1" ]
then
    sample="${SamplePath}/${SampleName}.csfasta"
    # Filtering reads that are not abundant enough for assembly
    if [ ! -f "${resultDir}/${SampleName}_filtered.csfasta" ] && [ -f "$sample" ] && [ "$abundance_filtering" -eq "1" ]
    then
        say "Filter reads depending on abundance"
        start_time=$(timer)
        fasta_seq="${resultDir}/${SampleName}.fasta"
        $csfastatofasta -i $sample -o $fasta_seq
        nb_reads=$(grep "^>" -c $sample)
        let "size=50*$nb_reads"
        $jellyfish count -C -m $kmer -t $NbProc -s $size -c 3 -o ${resultDir}/count_jellyfish --lower-count $lowerCount $fasta_seq
        nb_jelly=$(find ${resultDir}/ -name "count_jellyfish_*" |wc -l)
        if [ "$nb_jelly" -gt "1" ]
        then
            $jellyfish merge $(ls $output/count_jellyfish*) -o $2/count_jellyfish_merged
            jelly_output=${resultDir}/count_jellyfish_merged
        else
            jelly_output=${resultDir}/count_jellyfish_0
        fi
        $filterReads -i $fasta_seq -o ${resultDir}/${SampleName}_filtered.list -j $jelly_output -t $occurence
        $grabcataloguesequence -i ${resultDir}/${SampleName}_filtered.list -d $sample -o ${resultDir}/${SampleName}_filtered.csfasta
        # Sample is filtered with abundance
        sample="${resultDir}/${SampleName}_filtered.csfasta"
        if [ -f $fasta_seq ]
        then
            rm $fasta_seq
        fi
        if [ -f "$jelly_output" ]
        then
            rm $jelly_output
        fi
        say "Elapsed time to filter reads depending on abundance : $(timer $start_time)"
    fi
    # Filtering reads from others types than bacteria/fungi
    let "essai=1";
    if [ ! -f "${resultDir}/${SampleName}_${#filterRef[@]}.csfasta" ]
    then
        for db in ${filterRef[@]}
        do
            let "num=$essai-1";
            if [ ! -f "${resultDir}/${SampleName}_${essai}.csfasta" ] &&  [ -f "${resultDir}/${SampleName}_${num}.csfasta"  ] || [ "$essai" -ne "1" ]
            then
                say "Filter reads in $db"
                start_time=$(timer)
                # Next mapping
                $bowtie -f -C -S -v $NbMismatchMapping -p $NbProc -3 $readTrimLength $db ${resultDir}/${SampleName}_${num}.csfasta /dev/null --un ${resultDir}/${SampleName}_${essai}.csfasta -t  > ${logDir}/log_mapping_${SampleName}_${essai}.txt 2>&1
                check_file ${resultDir}/${SampleName}_${essai}.csfasta
                # Remove old file
                rm ${resultDir}/${SampleName}_${num}.csfasta
                say "Elapsed time to filter reads in $db : $(timer $start_time)"
            elif [ -f "$sample" ]  && [ "$essai" -eq "1" ] && [ ! -f "${resultDir}/${SampleName}_${essai}.csfasta"  ]
            then
                say "Filter reads in $db"
                start_time=$(timer)
                # First mapping
                $bowtie -f -C -S -v $NbMismatchMapping -p $NbProc -3 $readTrimLength $db $sample /dev/null --un ${resultDir}/${SampleName}_${essai}.csfasta -t  > ${logDir}/log_mapping_${SampleName}_${essai}.txt 2>&1
                check_file ${resultDir}/${SampleName}_${essai}.csfasta
                say "Elapsed time to filter reads in $db : $(timer $start_time)"
            fi
            let "essai=$essai+1";
        done
    fi

    # Triming
    filteredSample="${resultDir}/${SampleName}_${#filterRef[@]}.csfasta"
    if [ "$readTrimLength" -ne "0" ]
    then
        filteredSample="${resultDir}/${SampleName}_${#filterRef[@]}_trim_${readTrimLength}.csfasta"
    fi
    if [ -f  "${resultDir}/${SampleName}_${#filterRef[@]}.csfasta" ] && [ ! -f $filteredSample ]
    then
      say "Triming reads at $readTrimLength"
      start_time=$(timer)
      python $trimReads -f ${resultDir}/${SampleName}_${#filterRef[@]}.csfasta -s $readTrimLength -o $filteredSample 2> ${errorlogDir}/error_log_triming_${SampleName}.txt
      check_log ${errorlogDir}/error_log_triming_${SampleName}.txt
      check_file $filteredSample
      say "Elapsed time to trim : $(timer $start_time)"
    fi

    # Read assembly with denovo2
    if [ -f  "$filteredSample" ] && [ ! -f "${resultDir}/${SampleName}_nt_contigs.fa" ]
    then
        say "Assembly metagenome"
        start_time=$(timer)
        $denovo2/assemble.pl $filteredSample none $refLength -outdir ${resultDir}/${SampleName}_assembly -numcores $NbProc -cov_cutoff 3 -min_contig_lgth $contigLengthThreshold $denovo2Options > ${logDir}/log_assembly_${SampleName}.txt
        mv ${resultDir}/${SampleName}_assembly/nt_contigs.fa ${resultDir}/${SampleName}_nt_contigs.fa
        mv ${resultDir}/${SampleName}_assembly/analysis/nt_contigs/n50.stats.txt ${resultDir}/${SampleName}_n50_stats.txt
        mv ${resultDir}/${SampleName}_assembly/velvet/Log ${logDir}/log_velvet_${SampleName}.txt
        rm -r ${resultDir}/${SampleName}_assembly
        check_file ${resultDir}/${SampleName}_nt_contigs.fa
        say "Elapsed time to assembly : $(timer $start_time)"
    fi

    # Check if length filtering has to be done
    if [ "$contigLengthThreshold" -eq "300" ]
    then
        contigs="${resultDir}/${SampleName}_nt_contigs.fa"
    else
        contigs=${resultDir}/${SampleName}_nt_contigs.$contigLengthThreshold.fa
        # Filter the contigs of small size
        if [ -f "${resultDir}/${SampleName}_nt_contigs.fa" ] && [ ! -f "$contigs" ]
        then
            say "Filtering contigs of small size"
            start_time=$(timer)
            python $FilterFasta -f ${resultDir}/${SampleName}_nt_contigs.fa -t $contigLengthThreshold -s ${SampleName}  -o ${resultDir}/${SampleName}_nt_contigs.$contigLengthThreshold.fa 2> ${errorlogDir}/error_log_filterfasta_${SampleName}.txt
            check_log ${errorlogDir}/error_log_filterfasta_${SampleName}.txt
            check_file $contigs
            say "Elapsed time to filter contigs of small size : $(timer $start_time)"
        fi
    fi

     # Build bowtie on the assembled reads
    if [ -f "$contigs" ] && [ ! -f "$contigs.1.ebwt" ]
    then
        say "Bowtie-Build on the contigs"
        start_time=$(timer)
        $bowtie_build -f -C $contigs $contigs > ${logDir}/log_bowtiebuild_${SampleName}.txt 2> ${errorlogDir}/error_log_bowtiebuild_${SampleName}.txt
        check_log ${errorlogDir}/error_log_bowtiebuild_${SampleName}.txt
        say "Elapsed time to build bowtie reference on contigs : $(timer $start_time)"
    fi

    # Get unmapped reads
    if [ -f "$contigs.1.ebwt" ] && [ -f "${resultDir}/${SampleName}_${#filterRef[@]}.csfasta" ] && [ ! -f "${resultDir}/${SampleName}_unmapped.csfasta" ]
    then
        say "Get unmapped reads"
        start_time=$(timer)
        $bowtie -f -C -S -v $NbMismatchMapping -p $NbProc -3 $readTrimLength -t $contigs ${resultDir}/${SampleName}_${#filterRef[@]}.csfasta /dev/null --un ${resultDir}/${SampleName}_unmapped.csfasta -t > ${logDir}/log_unmapped_reads_${SampleName}.txt 2>&1
        check_file ${resultDir}/${SampleName}_unmapped.csfasta
        say "Elapsed time to get unmapped reads : $(timer $start_time)"
    fi

    # Predict genes
    if [ -f  "$contigs" ] && [ ! -f "${resultDir}/${SampleName}_gene.fna" ] && [ ! -f "${resultDir}/${SampleName}_prot.faa" ] && [ "$metagenemark" -eq "1" ]
    then
        say "Predict genes using Metagenemark"
        start_time=$(timer)
        $gmhmmp -a -d  -m $metagenemarkmod $contigs  -o ${resultDir}/${SampleName}.metagenemark 2> ${logDir}/log_metagenemark_${SampleName}.txt
        $extractMetagenemark -i ${resultDir}/${SampleName}.metagenemark -a ${resultDir}/${SampleName}_prot.faa -d ${resultDir}/${SampleName}_gene.fna
        check_file ${resultDir}/${SampleName}_gene.fna
        check_file ${resultDir}/${SampleName}_prot.faa
        say "Elapsed time to predict genes using Metagenemark : $(timer $start_time)"
    elif [ -f  "$contigs" ] && [ ! -f "${resultDir}/${SampleName}_gene.fna" ] && [ ! -f "${resultDir}/${SampleName}_prot.faa" ]
    then
        say "Predict genes using Prodigal"
        start_time=$(timer)
        $prodigal -i $contigs -a ${resultDir}/${SampleName}_prot.faa -d ${resultDir}/${SampleName}_gene.fna -p meta > ${resultDir}/${SampleName}.prodigal 2> ${logDir}/log_prodigal_${SampleName}.txt
        check_file ${resultDir}/${SampleName}_gene.fna
        check_file ${resultDir}/${SampleName}_prot.faa
        say "Elapsed time to predict genes using Prodigal : $(timer $start_time)"
    fi

    # Filter gene over their length
    if [ -f "${resultDir}/${SampleName}_gene.fna" ] && [ ! -f "${resultDir}/${SampleName}_gene_${geneLengthThreshold}.fna" ]
    then
        say "Remove genes shorter than $geneLengthThreshold nt"
        start_time=$(timer)
        python $FilterFasta -f ${resultDir}/${SampleName}_gene.fna -t $geneLengthThreshold -o ${resultDir}/${SampleName}_gene_${geneLengthThreshold}.fna 2> ${errorlogDir}/error_log_filterfasta_gene_${SampleName}.txt
        check_log ${errorlogDir}/error_log_filterfasta_gene_${SampleName}.txt
        check_file ${resultDir}/${SampleName}_gene_${geneLengthThreshold}.fna
        say "Elapsed time to remove short genes : $(timer $start_time)"
    fi

    # Remove genes that lack both ends
    if [ -f "${resultDir}/${SampleName}_gene_$geneLengthThreshold.fna" ] && [ ! -f "$input_gene" ]
    then
        say "Remove genes that lack both ends"
        start_time=$(timer)
        python $filterProdigal -i ${resultDir}/${SampleName}_gene_$geneLengthThreshold.fna -o $input_gene -c $filterGene 2> ${errorlogDir}/error_log_filterprodial_gene_${SampleName}.txt
        check_log ${errorlogDir}/error_log_filterprodial_gene_${SampleName}.txt
        check_file $input_gene
        # Extract list of element of interest
        grep ">" $input_gene |cut -d" " -f1 |cut -d">" -f2  >  ${resultDir}/${SampleName}_interest_list.txt
        say "Elapsed time to remove genes that lack both ends : $(timer $start_time)"
    fi

    # Filter proteins based on selected genes
    if [ -f "$input_gene" ] && [ -f "${resultDir}/${SampleName}_prot.faa" ] && [ -f "${resultDir}/${SampleName}_interest_list.txt" ] && [ ! -f "$input_protein" ]
    then
        say "Select proteins for which the gene has been selected"
        start_time=$(timer)
        python $extractProteins -i ${resultDir}/${SampleName}_interest_list.txt -d ${resultDir}/${SampleName}_prot.faa -o $input_protein 2> ${errorlogDir}/error_log_ExtractProteins_${SampleName}.txt  > ${logDir}/log_grab_sequence_${SampleName}.txt
        check_log ${errorlogDir}/error_log_ExtractProteins_${SampleName}.txt
        check_file $input_protein
        say "Elapsed time to select proteins for which the gene has been selected : $(timer $start_time)"
    fi
fi

if [ "$tax_annotation" -eq "1" ] || [ "$all" -eq "1" ]
then
    # Check for essential genes : MTU
    if [ -f $input_protein ] && [ ! -f ${resultDir}/${SampleName}_marker_cog.txt ]
    then
        say "Check for essential genes : MOTU"
        start_time=$(timer)
        mkdir ${resultDir}/fetch_result/
        $fetchmg -m extraction $input_protein -o ${resultDir}/fetch_result/ -d $input_gene -t $NbProc > ${logDir}/log_extract_marker_${SampleName}.txt
        wc -l ${resultDir}/fetch_result/temp/*.txt > ${resultDir}/${SampleName}_marker_cog.txt
        cat ${resultDir}/fetch_result/*.fna > ${resultDir}/${SampleName}_gene_marker.fna
        say "Elapsed time for MOTU : $(timer $start_time)"
    fi

    # Taxonomic annotation of genes with blastn on nt
    if [ -f "$input_gene" ] && [ ! -f "${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_nt.txt" ]
    then
        say "Taxonomic annotation with blastn on nt"
        start_time=$(timer)
        # -use_index true -db_soft_mask 11
        $blastn -query $input_gene -db $ntBlast -outfmt "6 qseqid sseqid qlen length pident evalue" -evalue $evalueTaxAnnot -num_threads $NbProc -out ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_nt.txt -max_target_seqs $maxTargetSeqs -task megablast 2> ${errorlogDir}/error_log_blastn_nt_${SampleName}.txt
        check_log ${errorlogDir}/error_log_blastn_nt_${SampleName}.txt
        check_file ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_nt.txt
        say "Elapsed time for taxonomic annotation with blastn on nt : $(timer $start_time)"
    fi

    # Extract annotations on ncbi
    if [ -f "${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_nt.txt" ] && [ ! -f "${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80.txt" ] && [ -f "$input_gene" ]
    then
        say "Analyze results on ncbi"
        start_time=$(timer)
        # ncbi annotation
        ncbi_annotation="${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_nt.txt"
        $gettaxonomy -i $ncbi_annotation -t $gitaxidnucl -n $names -d $nodes -o ${resultDir}/${SampleName}_taxonomy.txt 2> ${errorlogDir}/error_log_gettaxonomy_${SampleName}.txt
        check_file ${resultDir}/${SampleName}_taxonomy.txt
        check_log ${errorlogDir}/error_log_gettaxonomy_${SampleName}.txt
        python $extractNCBIDB -f $ncbi_annotation -g ${resultDir}/${SampleName}_taxonomy.txt -o ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80.txt -nb $numberBestannotation  -fc 80 2> ${errorlogDir}/error_log_extractncbidb_${SampleName}.txt
        check_log ${errorlogDir}/error_log_extractncbidb_${SampleName}.txt
        check_file ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80.txt
        say "Elapsed time to analyze results on ncbi : $(timer $start_time)"
    fi
fi


## Annotation
if [ "$func_annotation" -eq "1" ] || [ "$all" -eq "1" ]
then
    # blastp on EggNOG
    if [ -f $input_protein ] && [ ! -f ${resultDir}/${SampleName}_prot_eggnog_annotation.txt ]
    then
        say "Functional annotation with EggNOG V3"
        start_time=$(timer)
        $blastp -query $input_protein -db $eggnogRefBlast -out ${resultDir}/${SampleName}_prot_eggnog_annotation.txt -outfmt "6 qseqid sseqid qlen length pident evalue" -evalue $evalueFuncAnnot -num_threads $NbProc -max_target_seqs $maxTargetSeqs 2> ${errorlogDir}/error_log_blastp_eggnog_${SampleName}.txt
        check_log ${errorlogDir}/error_log_blastp_eggnog_${SampleName}.txt
        #check_file  ${resultDir}/${SampleName}_prot_eggnog_annotation.txt
        say "Elapsed time to get functional annotation with EggNOG V3 : $(timer $start_time)"
    fi

    # Extract eggnog annotations
    if [ -f  ${resultDir}/${SampleName}_prot_eggnog_annotation.txt ] && [ ! -f ${resultDir}/${SampleName}_prot_eggnog_annotation_filtered.txt ] && [ -f $input_protein ]
    then
        say "Extract EggNOG V3 result"
        start_time=$(timer)
        # -c -p $(dirname $clustalo)
        python $extractEggNog -b ${resultDir}/${SampleName}_prot_eggnog_annotation.txt -m $eggnogmemRef -n $eggnogdescRef -o ${resultDir}/${SampleName}_prot_eggnog_annotation_filtered.txt -f $eggnogfunccatRef -a $eggnogfunccats -nb $numberBestannotation 2> ${errorlogDir}/error_log_extractegg_${SampleName}.txt
        check_log ${errorlogDir}/error_log_extractegg_${SampleName}.txt
        check_file ${resultDir}/${SampleName}_prot_eggnog_annotation_filtered.txt
        say "Elapsed time to extract EggNOG V3 result : $(timer $start_time)"
    fi

    # hmmer on cazy
    if [ -f $input_protein ] && [ ! -f ${resultDir}/${SampleName}_prot_cazy_annotation.txt ]
    then
        say "Functional annotation with cazy"
        start_time=$(timer)
        $hmmscan --cpu $NbProc --tblout ${resultDir}/${SampleName}_prot_cazy_annotation.txt $cazyRefHmmer $input_protein > ${logDir}/log_cazy_hmm_${SampleName}.txt 2> ${errorlogDir}/error_log_cazy_hmm_${SampleName}.txt
        check_log ${errorlogDir}/error_log_cazy_hmm_${SampleName}.txt
        #check_file ${resultDir}/${SampleName}_prot_cazy_annotation.txt
        say "Elapsed time to get functional annotation with cazy : $(timer $start_time)"
    fi

    # Extract cazy annotations
    if [ -f ${resultDir}/${SampleName}_prot_cazy_annotation.txt ] && [ ! -f ${resultDir}/${SampleName}_prot_cazy_annotation_filtered.txt ]
    then
        say "Extract cazy result"
        start_time=$(timer)
        python $extractCazy -i ${resultDir}/${SampleName}_prot_cazy_annotation.txt -f $famInfoRef -o ${resultDir}/${SampleName}_prot_cazy_annotation_filtered.txt 2> ${errorlogDir}/error_log_extractcazy_${SampleName}.txt
        #check_log ${errorlogDir}/error_log_extractcazy_${SampleName}.txt
        #check_file ${resultDir}/${SampleName}_prot_cazy_annotation_filtered.txt
        say "Elapsed time to extract cazy result : $(timer $start_time)"
    fi
fi
say "Metagenome assembly is done. Elapsed time: $(timer $wall_time)"
