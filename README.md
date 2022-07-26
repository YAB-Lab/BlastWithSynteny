# BlastWithSynteny
Performs tBLASTn of a focal protein and its flanking neighbors against a set of target genomes. Return a GTF/GFF of the target hits and their corresponding query BLAST overlaps.


	Usage:  BlastSynteny [-g <GFF>] [-p <protein FASTA>]
		[-t <genomes/gff folder>] [-a <annotation_file>]
		[-i <gene_id>] ... [OPTIONS]


	-g <GFF>			The query GFF annotation file.
	-p <protein FASTA>		The query protein sequence database in FASTA format
	-f <genomes/gff folder>		The folder containing the target genomes and their GTF annotations.
					[file names should have extension '.fa' and '.gff']
	-a <annotation_file>		A tab separated file with three columns: gene_id, transcript_id, and protein_id
	-i <gene_id>			The query gene id as it appears in 'annotation_file'


	Optional arguments:
	-d <flank_dist>			The flanking region size in bp [default = 25000]


		Requirements: The following programs/scripts need to be in your path variable:

		1. fastagrep.pl
		2. FastaNamesSizes_syn.pl
		3. doSyntenyFilter.pl (needs to be installed as part of the 'SynBlast' perl package)
		4. gffread (obtain from 'https://github.com/gpertea/gffread', NOT conda!)
		5. NCBI BLAST+

		Other dependencies can be loaded from the 'synblast_env.yml' conda environment


