#!/bin/bash


## help menu
if [ "$1" == "-h" ] || [ $# -eq 0 ]; then
  echo -e "\n\t"
  pyfiglet -f puffy BlastSynteny
  echo -e "\n\tPerforms synteny BLAST between a focal query protein and its neighbours\n"
  echo -e "\n\t\tUsage: `basename $0` <query GFF> <query protein FASTA> <target genomes/gff folder> <query gene_id> <annotation_file>"
  echo -e "\n\n\t\t<query GFF>\t\t\tThe query GFF annotation file."
  echo -e "\t\t<query protein FASTA>\t\tThe query protein sequence database in FASTA format"
  echo -e "\t\t<target genomes/gff folder>\tThe folder containing the target genomes and their GTF annotations.\n\t\t\t\t\t\t[file names should have extension '.fa' and '.gff']"
  echo -e "\t\t<query gene_id>\t\t\tThe query gene id as it appears in 'annotation_file'"
  echo -e "\t\t<Annotation file>\t\tA tab separated file with three columns: gene_id, transcript_id, and protein_id\n\n"
  echo -e "\t\tRequirements: The following programs/scripts need to be in your path variable:\n"
  echo -e "\t\t\t1. fastagrep.pl"
  echo -e "\t\t\t2. FastaNamesSizes_syn.pl"
  echo -e "\t\t\t3. doSyntenyFilter.pl (needs to be installed as part of the 'SynBlast' perl package)"
  echo -e "\t\t\t4. gffread (obtain from 'https://github.com/gpertea/gffread', NOT conda!)"
  echo -e "\t\t\t5. BLAST+\n"
  echo -e "\t\t\tOther dependencies can be loaded from the 'synblast' conda environment\n\n"

  exit 0
fi

pyfiglet -f puffy BlastSynteny

# echo -e "\n\n\t Starting BlastSynteny... "
## make the temporary working directory
mkdir -p .tmp_synblast_folder.$4

pyfiglet -f digital Starting BlastSynteny...

echo -e "\n\n\t Obtaining the transcript sequences for gene $4"

num_transcrits=$(grep $4 $5 | cut -f 2 | wc -l)

grep $4 $5 | cut -f 2 -> .tmp_synblast_folder.$4/transcript_list

echo -e "\n\n\t $4 has $num_transcrits transcript(s):\n"

cat .tmp_synblast_folder.$4/transcript_list | sed 's/^/\t\t/g'

transcript=$(head -1 .tmp_synblast_folder.$4/transcript_list)

echo -e "\n\n\t\t Running the analysis on the first transcrtipt: $transcript"

## Starting process. Define GTF elements
query_chr=$(grep $transcript $1 | awk '$3 == "transcript"' | cut -f 1)
query_min=$(grep $transcript $1 | awk '$3 == "transcript"' | cut -f 4)
query_max=$(grep $transcript $1 | awk '$3 == "transcript"' | cut -f 5)

echo -e "\n\n\t\t Query gene coordinates: $query_chr $query_min $query_max"

## create the temporary protein info file
awk -v chr="$query_chr" '$1 == chr' $1 | \
  awk -v min="$query_min" '$4 > min-25000' | \
  awk -v max="$query_max" '$4 < max+25000' | \
  awk '$3 == "transcript"' | \
  cut -f 4,5,7,9 | \
  grep protein_coding | \
  awk '{print $1, $2, $7, $5, $9, $3}' | \
  sed 's/"//g' | \
  sed 's/;//g' | \
  sed 's/+$/1/g' | \
  sed 's/-$/-1/g' | \
  tr ' ' '\t' > .tmp_synblast_folder.$4/tmp_info_file.txt

## extract the protein sequences of the focal gene and its neighbours
cut -f 3 .tmp_synblast_folder.$4/tmp_info_file.txt | \
  fastagrep.pl -f - $2 > .tmp_synblast_folder.$4/$4.proteins.fasta


echo -e "\n\n\t\t# The query FASTA file contains:"
## output the protein lengths:
FastaNamesSizes_syn.pl .tmp_synblast_folder.$4/$4.proteins.fasta > .tmp_synblast_folder.$4/$4.proteins.lengths

## generate the complete gene/transcript/protein IDs file:
grep ">" .tmp_synblast_folder.$4/$4.proteins.fasta | \
  sed 's/pep primary.*gene:/\t/g' | \
  sed 's/ transcript:/\t/g' | \
  sed 's/ gene_biotype.* gene_symbol:/\t/g' | \
  sed 's/ description.*//g' | \
  sed 's/>//g' > .tmp_synblast_folder.$4/$4.annot

## combining lengths with protein info file
join .tmp_synblast_folder.$4/$4.proteins.lengths .tmp_synblast_folder.$4/$4.annot | awk '{print $4, $1, $2, $3, $5}' | tr ' ' '\t' > .tmp_synblast_folder.$4/tmpFile

## Missing a couple of variables for the header:
strand_ori=$(grep $transcript .tmp_synblast_folder.$4/tmp_info_file.txt | cut -f 6)
prot_length=$(grep $transcript .tmp_synblast_folder.$4/tmpFile | cut -f 3)

## Create the header for the protein info file
echo "#protein.info file
#!organism=Drosophila_melanogaster
#!ensemblRelease=dec2006
#!coordsystem=chromosome
#!seqregion=CHR
#!from=PROT_START
#!to=PROT_END
#!strand=STRAND_ORI
#!length=PROT_LENGTH
#!focalgene=GENE_NAME
#!flankingsize=5e5
#!columns=startPos,endPos,meanPos,protID,protName,ori,geneID,geneName,transID,protLength" | \
  sed "s/CHR/$query_chr/g" | \
  sed "s/PROT_START/$query_min/g" | \
  sed "s/PROT_END/$query_max/g" | \
  sed "s/STRAND_ORI/$strand_ori/g" | \
  sed "s/PROT_LENGTH/$prot_length/g" | \
  sed "s/GENE_NAME/$4/g" > .tmp_synblast_folder.$4/$4.header

## Make the final protein info file:
cat .tmp_synblast_folder.$4/tmp_info_file.txt | \
    awk '{print $3, $1, $2, $4, $5, $6}' | \
    join - .tmp_synblast_folder.$4/tmpFile | \
    awk 'BEGIN{OFS="\t"}{print $2, $3, ($2 + $3)/2, $7, $5, $6, $4, $5, $1, $8}' | \
    cat .tmp_synblast_folder.$4/$4.header - > .tmp_synblast_folder.$4/$4.proteins.info

###
###

## prep the reference genome BLAST databses:
mkdir .tmp_synblast_folder.$4/blastdbs
cp $3/*.fa .tmp_synblast_folder.$4/blastdbs

# make a list of species
ls -1 $3 | grep ".fa$" | sed 's/\..*//g' > .tmp_synblast_folder.$4/ref.list

##
echo -e "\n\n\t Building reference BLAST databases:"

for genome in `cat .tmp_synblast_folder.$4/ref.list`
  do
    echo -e "\n\t Building reference for $genome genome"
    makeblastdb -in .tmp_synblast_folder.$4/blastdbs/$genome.fa -dbtype nucl >> .tmp_synblast_folder.$4/blastDB.log 2>&1
  done

## Run blast processes for each reference genome
mkdir .tmp_synblast_folder.$4/tblastn/

for genome in `cat .tmp_synblast_folder.$4/ref.list`
  do
    echo -e "\n\t Running tBLASTn against $genome"
    tblastn \
    -query .tmp_synblast_folder.$4/$4.proteins.fasta \
    -db .tmp_synblast_folder.$4/blastdbs/$genome.fa \
    -outfmt 6 \
    -evalue 1e-5 > .tmp_synblast_folder.$4/tblastn/$4.proteins.info.$genome.blastresult
  done

## Now run synteny filter

echo -e "\n\n\n"

pyfiglet -f digital Filtering synteny hits...

for genome in `cat .tmp_synblast_folder.$4/ref.list`
  do
    doSyntenyFilter.pl \
    -d$genome \
    -r.tmp_synblast_folder.$4/tblastn/ \
    -q40 \
    -f4 \
    -a$4.proteins.info \
    -o \
    -g10 \
    -n10 \
    .tmp_synblast_folder.$4/$4.proteins.info
  done


# doSyntenyFilter.pl -d$genome -r.tmp_synblast_folder.$4/tblastn/ -q40 -f4 -a$4.proteins.info -o -g10 -n10 .tmp_synblast_folder.$4/$4.proteins.info


echo -e "\n\n\n"

pyfiglet -f digital Creating GTF outputs

## Create a GTF file of the reference hits
for genome in `cat .tmp_synblast_folder.$4/ref.list`
  do
    sed '/^#/d' .tmp_synblast_folder.$4/$4.proteins.syntenyRegions/$genome.*.blastresult | \
    awk -v OFS='\t' 'NR == FNR{a1[$1]=$2; a2[$1]=$3; a3[$1]=$4; next}; $1 in a1{$13=a1[$1]; $14=a2[$1]; $15=a3[$1]};{print}' .tmp_synblast_folder.$4/$4.annot - | \
    awk 'BEGIN{FS=OFS="\t"}{print $2, "synBlst", "BLST", $9, $10, ".", ".", ".", "transcript_id \""$14"\"; gene_id \""$13"\"; gene_name \""$15"\"; perc_ident \""$3"\"; evalue \""$11"\""}' > .tmp_synblast_folder.$4/$4.$genome.blast_results.gtf
  done


## Create part of GTF for BLAST hits region:
# get chromosome ID:
for genome in `cat .tmp_synblast_folder.$4/ref.list`
  do
    sed '/^#/d' .tmp_synblast_folder.$4/$4.proteins.syntenyRegions/$genome.*.blastresult | \
    cut -f 2 |\
    sort -u > .tmp_synblast_folder.$4/$4.$genome.chr

    sed '/^#/d' .tmp_synblast_folder.$4/$4.proteins.syntenyRegions/$genome.*.blastresult | \
    cut -f 9 |\
    sort -n |\
    head -1 > .tmp_synblast_folder.$4/$4.$genome.min

    sed '/^#/d' .tmp_synblast_folder.$4/$4.proteins.syntenyRegions/$genome.*.blastresult | \
    cut -f 10 |\
    sort -n |\
    tail -1 > .tmp_synblast_folder.$4/$4.$genome.max

    cat .tmp_synblast_folder.$4/$4.$genome.chr \
        .tmp_synblast_folder.$4/$4.$genome.min \
        .tmp_synblast_folder.$4/$4.$genome.max | \
        tr '\n' '\t' | \
        sed 's/$/\n/g' > .tmp_synblast_folder.$4/$4.$genome.range

    rm .tmp_synblast_folder.$4/$4.$genome.chr \
        .tmp_synblast_folder.$4/$4.$genome.min \
        .tmp_synblast_folder.$4/$4.$genome.max
  done



for genome in `cat .tmp_synblast_folder.$4/ref.list`
  do
    desired_range=$(cat .tmp_synblast_folder.$4/$4.$genome.range | \
      sed 's/\t/:/1' | \
      sed 's/\t/../1')
    gffread references/$genome.gff \
      -r $desired_range \
      -T |\
      cat - .tmp_synblast_folder.$4/$4.$genome.blast_results.gtf > .tmp_synblast_folder.$4/$genome.$4.SynRegion.gtf
  done


### Pass final outputs
mkdir -p $4.synteny_blast_results

for genome in `cat .tmp_synblast_folder.$4/ref.list`
  do
    cp .tmp_synblast_folder.$4/$genome.$4.SynRegion.gtf $4.synteny_blast_results
  done

pyfiglet -f digital FIN!

## clean up
# rm -r .tmp_synblast_folder.$4
