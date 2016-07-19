# ------------------------------------------------------------------
# SOMA version 0.0
# Authors: Amine Ghozlane (amine.ghozlane@jouy.inra.fr)
#          Mathieu Almeida
#          Florian Plaza Oñate
# Special thanks : Anne-Sophie Alvarez
# Title:  Assemble_SOLiD_lite-metagenome
# Description : Assembly pipeline for SOLID reads with annotation
# ------------------------------------------------------------------

#########################
# Activate Metagenemark #
#########################

# Download Metagenemark 64 bits version :
http://exon.gatech.edu/genemark/license_download.cgi

# Copy mgm directory (inside the archive) in soma/MetaGeneMark/

# Activate MetaGeneMark
mv ./soma/MetaGeneMark/mgm/gm_key  $HOME/.gm_key


#################
# Data for SOMA #
#################

# Download and move data.tar.gz in soma directory, then :
```
wget http://mgps.eu/ibs-tools/data.tar.gz
tar -zxf data.tar.gz
mv data soma/
cd soma/data/
./install_databank.sh
```

############
# Run SOMA #
############

# Help
```
./soma/soma.sh
Use '--help' to print detailed descriptions of command line arguments
```

# Usage :
# case 1 - Assembly: 
`./soma.sh -i <read.csfasta> -o </path/to/result/directory/> --assembly (--hsize 23 -r 20000000 for 35nt reads or --hsize 25 -r 100000000 for 50nt reads --abundance_filtering --lowerCount 3 --occurence 3 --kmer 15  for low quality sequencing)`

# case 2 - Taxonomic annotation/MOTU: 
`./soma.sh -g <gene.fasta> -o </path/to/result/directory/> --tax_annotation`

# case 3 - Functional annotation:
`./soma.sh -p <protein.fasta> -o </path/to/result/directory/> --func_annotation`

# all :
`./soma.sh -i <read.csfasta> -q <rea.quality> -o </path/to/result/directory/>`

##############################
# Particular case : Taxonomy #
##############################

Python package generated for get_taxonomy might not fit to every system
configuration.
To test it, try :

`./soma/get_taxonomy/get_taxonomy -h`

# Result :
```
usage: ./soma/get_taxonomy -h (see also ftp://ftp.ncbi.nih.gov/pub/taxonomy/)

Get NCBI taxonomy with lineage.

optional arguments:
  -h, --help            show this help message and exit
  -i BLAST_OUTPUT_FILE  blast_output_file
  -t NCBI_TAXID_FILE    ncbi_taxid_file
  -n NCBI_NAMES_FILE    ncbi_names_file
  -d NCBI_NODES_FILE    ncbi_nodes_file
  -l LINEAGE            List of ranks to reported in taxonomy lineage (default
                        ['superkingdom','kingdom',
                        'phylum','class','order','family','genus','species'])
  -o TAXONOMY_FILE      Output taxonomy_file
```

# If failed, run as root :
`pip install cogent`
# In soma.sh, make the following changes (line 224-225) :
`gettaxonomy="$SCRIPTPATH/get_taxonomy/get_taxonomy.py"`

# Example dataset
One example dataset is available [here](http://mgps.eu/ibs-tools/example_soma.tar.gz)
