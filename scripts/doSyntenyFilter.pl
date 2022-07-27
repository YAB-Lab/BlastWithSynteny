#!/usr/bin/env perl
use SynBlast;
our $PROGN=[split(/\//, $0)]->[-1];
$VERSION = "SynBlast-" . $VERSION . "\n${PROGN}'s Time-stamp: <2008-08-19 14:59:27 joe>\nsee ${url} for updates.\n";
#######################################################
#doSyntenyFilter.pl  - part of the SynBlast package for Assisting the Analysis of Conserved Synteny Information
#Copyright (C) 2008  Joerg Lehmann
# 2006-01-29
# joe@bioinf.uni-leipzig.de
#######################################################

=head1 NAME

doSyntenyFilter.pl - a tool to filter/separate blast hits into regions of possibly conserved synteny compared to a reference set of sequences (query proteins).
Query proteins' annotation must be provided in a reference annotation file in tab-separated format (.protein.info -file, see bottom).
The output is a set of  tabular output files each containing only blast hits within a restricted target sequence range as specified.
Those sets of syntenic blast hits are additionally restricted by various parameters,
most importantly: each result set must contain at minimum a user-specified number of query-specific blast hits (-q/-Q options).

doSyntenyFilter.pl is part of the L<SynBlast> package for Assisting the Analysis of Conserved Synteny Information.

=head1 SYNOPSIS

B<doSyntenyFilter.pl>
[B<-d I<DBname>>]
[B<-e I<Ensembl_DB_Version_String>>]
[B<-r I<Blastresult_Source_Directory>>]
[B<-a I<BlastFilePraefix>>]
[B<-q I<MinQueryPercent>>]
[B<-Q I<MinQueryNumber>>]
[B<-f I<SizeFactor>>]
[B<-F I<MaxAbsoluteSize>>]
[B<-w I<Region_Destination_Directory>>]
[B<-o  (WriteResultFlag) >]
[B<-n I<MaxResultsPerCtgToWrite>>]
[B<-g I<MaxResultsToWrite>>]
[B<-c I<ScoreThreshold>>]
[B<-i I<MaxOverlapPercent>>]
[B<-W (WU-BLAST_Flag) >]
[B<-RI<chr/contig>>]
[B<-s  (disableUseOfScoreFlag) >]
[B<-G I<SlidingWindowStepSizeParam>>]

[B<--help>]
[B<--version>]
protein.info-File

=head1 REQUIRES

require MyUtils;
require MySyntenyUtils;
require Getopt::Std;
require Pod::Usage;

=head1 DESCRIPTION

The extraction of syntenic blast hit subsets w.r.t. a reference/query set is done via a score-guided sliding window approach,
aimed to find all maximal target region intervals (and their associated/contained blast hits) that fullfill the user-specified
constraints upon gene content (-q/-Q options) and maximal genomic extension (-f/-F options).
Further options described below.


=head1 OPTIONS

=over 5

=item B<-d> I<DBname> [String]

 To specify the name of the target species (or database name).
 This item is used to select the appropriate blast result files of type "[QueryID|<BlastFilePrefix>].<DBname>.blastresult".
 default = Homo_sapiens

=item B<-e> I<Ensembl_DB_Version_String> [String]

 To specify the Ensembl version of the target species database, e.g. "aug2006" for the Ensembl release 46 (archive string).
 This is optional but recommended if Ensembl databases are used that are not of the same release as the API version used (not recommended).
 Evaluated target regions are then correctly linked to the entries of the appropriate Ensembl archive web site.
 Note: Database versions of reference and target species data as well as the used Ensembl API version should usually be the same!
 The default is the archive string corresponding to the Ensembl API version in use.

=item B<-r> I<Blastresult_Source_Directory> [String]

 To specify the directory containing the appropriate blast result files, e.g. "/scratch/joe/tblastn"
 The format of file names for the blast result files has to be either "QueryID.<DBname>.blastresult" (by default)
 or "<BlastFilePraefix>.<DBname>.blastresult" (if <BlastFilePraefix> is set).
 default = "tblastn"

=item B<-a> I<BlastFilePraefix> [String]

 To specify an alternative blast file praefix and not to look for query specific blast result files,
 but a file named "<BlastFilePraefix>.<DBname>.blastresult"
 default = undef
 (There is one query specific blast result file for every single query sequence, named "QueryID.<DBname>.blastresult"

=item B<-q> I<MinQueryPercent> [Integer]

 To specify the Minimum-Query-Number constraint in percent (i.e. the number of query-specific HSPs to be contained
 in a reported candidate region at minimum, given in percent relative to the number of total query proteins in the reference set).
 Note that the maximum of both resulting values (L<MinQueryPercent> and L<MinQueryNumber>) is used as final QueryNumber constraint!
 default = 0
 (i.e. at least hits from L<MinQueryNumber> different queries are necessary for a region to be reported.)

=item B<-Q> I<MinQueryNumber> [Integer]

 To specify the Minimum-Query-Number constraint as absolute number (i.e. the number of query-specific HSPs to be contained
 in a reported candidate region at minimum).
 Note that the maximum of both resulting values (L<MinQueryPercent> and L<MinQueryNumber>) is used as final QueryNumber constraint!
 default = 2
 (i.e. hits from at least two different reference proteins are necessary for a region to be reported.)


=item B<-f> I<SizeFactor> [Float]

 To specify the maximal length allowed for a candidate syntenic target region, relative to the reference block's size.
 Depending on the size differences between reference and target genome and its distribution of genes,
 this parameter should be adapted.
 default = 2.0
 (i.e. twice of the original reference's size, depending on the positions/range of genes, is used to restrict
 the length of reported synteny regions)

=item B<-F> I<MaxAbsoluteSize> [Integer]

 To specify the maximal length allowed for a candidate syntenic target region in absolute numbers (in nts).
 Note that the minimum of both resulting sizes (L<SizeFactor>-based size and L<MaxAbsoluteSize>)
 is used as the final MaxLength constraint
 default = 0
 (i.e. not set/disabled)

=item B<-w> I<Region_Destination_Directory> [String]

 To specify where to save the extracted candidate regions (blast result files).
 default is the subdirectory with the name of the <protein.info-file> + ".syntenyRegions"

=item B<-o>

 Flag to enable the output of all candidate regions to files in L<Region_Destination_Directory>.
 default = not set
 (i.e. by default no output to any files)

=item B<-n> I<MaxResultsPerChromosomeToWrite> [Pos. Integer]

 To specify the number of regions for each target chromosome/contig to be reported at maximum.
 default = 100
 (i.e. only the first best 100 candidate regions per chromosome/contig, sorted according to the sum of
 maximum bitscores of all HSPs for each different query, are written to files/displayed)

=item B<-g> I<MaxResultsInTotalToWrite> [Pos. Integer]

 To specify the number of regions in total to be reported at maximum. Sorting of candidate regions is again
 according to the sum of query-specific bitscores (unless old version, with option -s is used).
 Note: This option has no function if the non-score version (option -s) is used.
 default = 0
 (i.e. no restriction, only the number of regions per chromosome/contig is restricted then)


=item B<-c> I<ScoreThreshold> [Integer]

 To specify the bit-score value that a single HSP must reach to be considered for the score-based selection
 e.g. 50. In general there is no need to change the default. This option could be used to eliminate regions with spurious hits,
 still experimental.
 default = 0

=item B<-i> I<MaxOverlapPercent> [Integer]

 To specify the maximally allowed overlap size in percent for lower-scoring overlapping result regions.
 Regions overlapping more than this fraction in size with an already/previously extracted candidate region
 are simply ignored.
 default = 25
 (i.e. up to 25 percent of reported candidate regions' sizes may be overlapping by default)

=item B<-W>

 Flag to use tabular blastresults in WU-BLAST format (different columns and orders).
 Note: output (synteny candidate regions) are always in a NCBI-BLAST-like tabular format!

=item B<-R> I<Chromosome/Contig> [String]

 To restrict the search/extraction on only a certain chromosome/contig [experimental]
 default = undef

=item B<-s>

 Flag to disable the use of score-based selection.
 (use of old filter version instead which does not take any score values into account; not recommended!)



=item B<-G> I<SlidingWindowStepSizeParam> [Integer]

 To set indirectly the stepSize for the sliding window method: referenceSize/param => stepSize,
 i.e. high values for this parameter increase the number of region-windows calculated for each target sequence
 default = 100. [experimental]


=head1 ARGUMENTS

=over 5

=item B<protein.info-file> [String]

 The reference info file to use.
 It specifies the annotation file containing the sequence features (proteins) to be used as reference syntenic region.
 For Ensembl genomes, L<getEnsemblProteins.pl> performs the retrieval/creation of such a file automatically.
 The format is a tab-separated list of features describing each protein-coding gene of the set, one per line.
 Definition of columns (column names separated by comma) and some more annotation descriptors have to be specified in the head,
 see the following example:

  # protein.info file manually generated
  #!organism=Branchiostoma
  #!coordsystem=scaffold
  #!seqregion=scaffold_89
  #!from=1
  #!to=650000
  #!strand=1
  #!focalgene=89000016
  #!columns=startPos,endPos,geneName,protID,ori,geneID,numExons,protLength
  42268   52505   gene1   fgenesh2_pm.scaffold_89000001   -1      89000001        3       592
  136667  137756  gene2   fgenesh2_pg.scaffold_89000006   1       89000006        2       121
  139176  141632  -       fgenesh2_pg.scaffold_89000007   -1      89000007        5       452
  147110  150418  -       fgenesh2_pg.scaffold_89000009   -1      89000009        1       1103
  #testcomment
  163582  175636  -       fgenesh2_pg.scaffold_89000011   1       89000011        15      649
  ...


=head1 AUTHOR

Joerg Lehmann, <joe@bioinf.uni-leipzig.de>

=head1 AVAILABILITY

http://www.bioinf.uni-leipzig.de/Software/SynBlast/


=head1 LICENSE

This  program is free  software; you can redistribute  it and/or
modify it  under the terms of the GNU  General Public License as
published by the  Free Software Foundation; either version 2, or
(at your option) any later version.

This program is  distributed in the hope that it will be useful,
but WITHOUT  ANY WARRANTY; without even the  implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You can view the  GNU General Public License, online, at the GNU
Project's homepage; see <http://www.gnu.org/licenses/gpl.html>.

See also L<perlgpl> for more information.

=head1  SEE ALSO

L<SynBlast::MySyntenyAln>,
L<SynBlast::MySyntenyUtils>,
L<SynBlast::MyUtils>,
L<getEnsemblProteins.pl>,
L<doBlastJobs.pl>
L<doSyntenyAlignment.pl>,
L<getBestBlastHit.pl>

=cut


use strict;
use warnings;

use SynBlast::MyUtils;
use SynBlast::MySyntenyUtils;
use Getopt::Std;
use Pod::Usage;
use DataBrowser;

$|=1;

############### options ###################################################################
my $cmdlineWas = ${main::PROGN} . " " . join(" ", @ARGV);

$Getopt::Std::STANDARD_HELP_VERSION = 1;
my %opts=();
getopts('d:a:q:Q:of:F:r:w:n:g:e:sc:i:WG:R:', \%opts) || pod2usage(-verbose => 0);
pod2usage(-verbose => 0, -message => "$PROGN: Error - exactly one InfoFile needed.\n") unless @ARGV==1;

my $infoFile=$ARGV[0];
my $infoName=substr($infoFile, 0, -5); # infoFile is assumed to be with ending .info
my $DBname="Homo_sapiens"; ## default database name is Homo_sapiens
$DBname = $opts{"d"} if(exists($opts{"d"}));

my $wublastFl=0;
$wublastFl = 1 if(exists($opts{"W"}));

my $restrChr=0;
$restrChr=$opts{"R"} if(exists($opts{"R"}));

my $EnsemblD=undef; ##"dec2006";  ## to connect the blast-hits with its Ensembl release string (manually to be selected)
$EnsemblD = lc($opts{"e"}) if(exists($opts{"e"}));
# $EnsemblD = lc(SynBlast::EnsEMBLaccess::getRelMonthString()) unless($EnsemblD);

my $maxPercOverlap=25;
$maxPercOverlap=int($opts{"i"}) if(exists($opts{"i"}));

my $writeFlag=0;# default behaviour: only show info, do not save regions as files
$writeFlag=int($opts{"o"}) if(exists($opts{"o"}));

my $scOrderFl=0;
$scOrderFl=int($opts{"s"}) if(exists($opts{"s"}));

my $mScThreshold=0;
$mScThreshold=int($opts{"c"}) if(exists($opts{"c"}));

my $stepGrain=100.;   ## the stepsize for the sliding window (on target-nt-basis) is then refsize/stepGrainSize (a higher step grain size results in smaller stepsizes, making the window more accurate w.r.t. scoreMaximaization
$stepGrain = int($opts{"G"}) if(exists($opts{"G"}));

my $gFl=0; # max number of regions to write out
if(exists($opts{"g"})) {
	$gFl=(int($opts{"g"}));
	die "Error: -g option needs positive value\n" unless ($gFl>0);
}

my $nFl=100; # max number of regions per chr to write out
$nFl=(int($opts{"n"})) if(exists($opts{"n"}));
die "Error: -n option needs positive value\n" unless ($nFl>0);

# default flanking size factor is 4 (i.e. 4x the effective length of the flanked sequence region around query protein is used as limit for size of regions in target
my $flP=2;
$flP=$opts{"f"} if(exists($opts{"f"}));
my $flSize=0;
$flSize=$opts{"F"} if(exists($opts{"F"}));

my $bSrcDir = "tblastn"; # default source dir (tblastn-result-files) is ./tblastn
my $DirBase= qx(pwd); chomp($DirBase);
$bSrcDir=$opts{"r"} if(exists($opts{"r"}));
SynBlast::MyUtils::checkRelDir(\$bSrcDir, $DirBase);

# default result directory for filtered blast hits (syntenyRegions) is "INFOFILEprefix.syntenyRegions"
my $synResultDir= $infoName . ".syntenyRegions";
$synResultDir=$opts{"w"} if(exists($opts{"w"}));
SynBlast::MyUtils::checkRelDir(\$synResultDir, $DirBase);

SynBlast::MyUtils::checkStrEnd(\$bSrcDir, "/"); # && print "blastresult directory name changed to " . $bSrcDir . "\n";
SynBlast::MyUtils::checkStrEnd(\$synResultDir, "/"); # && print "filter result directory name changed to " . $synResultDir . "\n";

if( -d "${bSrcDir}" ) {
    print "# using ${bSrcDir} as source directory for the blastresult files...\n";
} else {
    die "Error: directory \"${bSrcDir}\" (source for the blastresult_files, -r option) does not exist!\n";
}

if($writeFlag && (! (-d "${synResultDir}" ))) {
    mkdir( "${synResultDir}" , 0755) or die "Error: Could not create directory \"${synResultDir}\" (-w option, syntenyResult directory): $!\n";
}
print STDERR "# using ${synResultDir} as output directory for the synteny-filtered blasthits...\n" if($writeFlag);

## creating list of exons out of infoFile...
my @exons=();
my %exonIDIdxHash=();
my @refExons=();
my @exonList=();
my($icolIdx, $infoH, $realFlankS, $realSlength) = SynBlast::MySyntenyUtils::infofile2Exonlist("${infoFile}", \@exons, \%exonIDIdxHash, \@refExons, \@exonList, 0, "protID");
my $exonNumb = @exonList;
my @fileList= @exonList;
print STDERR "${infoFile} read. (" . @fileList . " entries)\n";

### one alternative blastfile to be used, instead of single blastfiles per QueryID
if(exists($opts{"a"})) {
    @fileList=( $opts{"a"} );
}

# minimum number of different Queries in synteny-FilterResult (in percent, but at least 2 different queries are default)
my $minQP=0;
$minQP=int($opts{"q"}) if(exists($opts{"q"}));
my $minQN=2;
$minQN=int($opts{"Q"}) if(exists($opts{"Q"}));
$minQN=SynBlast::MyUtils::max($minQN, int($exonNumb * $minQP / 100.0));
#print "# Query-Number constraint (QN) is set to ${minQN}.\n";

my %blastHits=();
my $scount=0;
my $totalFastas=@fileList;
foreach my $file (@fileList) {
    ++$scount;
    my $filestr = "${bSrcDir}/${file}.${DBname}.blastresult";
    print STDERR "filtering blastfile \"${file}\" (${scount}/${totalFastas})...\n";
    filterBlastFile2HashEq(\%blastHits, $filestr, \%exonIDIdxHash, undef, $wublastFl);
}

if($flSize) {
	$flSize = SynBlast::MyUtils::min($flSize, int($flP * $realSlength));
} else {
	$flSize = int($flP * $realSlength);
}

my $flSizeStr= SynBlast::MyUtils::formatMyNumber($flSize, "'", 3);
print "\n" . ("#" x40) . "\n"
			,"# Probable Syntenic Regions\n"
			,"# -of maximum size "
			, ( $flSize == int($flP * $realSlength)) ? ("${flP} x ${realSlength} nt (" . $flSizeStr . " nt)\n") : "${flSizeStr} nt\n"
			,"# -containing at least ${minQN} hits of different query\n"
			, ("#" x40) . "\n";

			#,"chromosome\tregionID\tsaveFlag\n############################\n";

if(! $scOrderFl) {
	print "Regions sorted by SumOfMaxBitScoresPerQuery (and secondarily by totalScore) (only up to those with fullfilled QN-constraint, skipping overlapping ones):\n";
	print "SumOfMaxScores\tstart\trstart\tend\tQN\ttotalHSPcount\ttotalscore\t%ofMaxRange\n";

	## start processing of chromosomes in order of total hits per chromosome
	my %blHitChrCounts=();
	foreach  my $cchr (keys %blastHits) {
		push(@{$blHitChrCounts{keys(%{$blastHits{$cchr}})}}, $cchr);
	}
	my %allRegionsBySc=();
	my $savedFiles="";
	foreach  my $cnt (sort {$a <=> $b} (keys %blHitChrCounts)) {
		if($cnt < $minQN) {
			print "." x scalar(@{$blHitChrCounts{$cnt}});
			next;
		}
		foreach  my $cchr (@{$blHitChrCounts{$cnt}}) {

			next if($restrChr && ($restrChr ne $cchr));


			#calculate SlidingWindow-value to identify regions
			my $distribution=calcSlidingWindowStat($blastHits{$cchr},  $flP,  $realSlength, $minQN, $mScThreshold, $flSize, $stepGrain);
			next unless(defined($distribution));



			my %regionsByMaxSc=();


			foreach my $ii (keys %{$distribution}) {



				my $percOfMaxRange=int(($distribution->{$ii}->[3] - $distribution->{$ii}->[6] +1)/($flSize) * 1000)/10.;

				### now use as second order-by field/criteria the total sum of scores (for regions of same max-score-sum-per-query)
				 ### 																* indicates where information is w.r.t. only those hits in a region,
				 ###																	that are above a scoreThreshold ($mScThreshold)
				 ### 					SumOf_MaxScorePerQuery*,   totalScoreSum*
			### 					diffQN   SumOf_MaxScorePerQuery*,   totalScoreSum*
				#push(@{$regionsByMaxSc{$distribution->{$ii}->[2]}->{$distribution->{$ii}->[4]}->{$distribution->{$ii}->[1]}},  [ (
				push(@{$regionsByMaxSc{"1"}->{$distribution->{$ii}->[4]}->{$distribution->{$ii}->[1]}},  [ (
									$ii  											#  genom. start (of sl.window)
									,  $distribution->{$ii}->[6]  						# genomic start used,
									,  $distribution->{$ii}->[3]						# genomic end used,
									, $distribution->{$ii}->[2]						#number of diff. queries *
									,  $distribution->{$ii}->[0]						#total number of HSPs contained
									, $distribution->{$ii}->[1]  						## total score sum *
									,  $percOfMaxRange							### percent of maxRangeUsed
									) ]);
			}
			#next;
			my $outstr="";
			my $breakFl=0;
			my @doneIntervs=();
			## sort order of evaluation is decreasing, first in SumOfMaxScoresPerQuery, second by totalScoreSum

			foreach my $QNkey (reverse (sort {$a <=> $b} (keys %regionsByMaxSc) ) ) {
				foreach my $scKey (reverse (sort {$a <=> $b} (keys %{$regionsByMaxSc{$QNkey}} ) ) ) {
					foreach my $tscKey (reverse (sort {$a <=> $b} (keys %{$regionsByMaxSc{$QNkey}->{$scKey}} ) ) ) {
						foreach my $entr  (@{$regionsByMaxSc{$QNkey}->{$scKey}->{$tscKey}}) {
							next if($entr->[3] < $minQN); ## QN not reached?

						### check on overlap with previous interval set
						my $oovlp=0;
						foreach my $intvl (@doneIntervs) {
							my $ooolp= (SynBlast::MyUtils::getOverlapSize(@{$intvl},  $entr->[1], $entr->[2]));
							$oovlp +=1 if($ooolp > int($maxPercOverlap * (abs($entr->[2]-$entr->[1]) + 1) / 100. ));
							### only skip overlapping candidates if overlapping more than maxPercOverlap% in length
						}
						next if($oovlp);	#print "/";

						my $regNaame=$scKey . "\t" . join("\t", @{$entr});
						#$outstr .= $regNaame;
						$outstr .= #$QNkey . "\t" .
								$regNaame;

						if((@doneIntervs < $nFl)) {  ### save region to syntenyRegion-File...
							if(! $gFl) {
								if($writeFlag) {
									$outstr .= "\t(saved)\n";
									# write region to file
									$savedFiles .= writeSynRegToFile(@doneIntervs +1,  $regNaame, $distribution->{$entr->[0]}->[5], $entr, $cchr,  $realSlength, $flP, undef, $flSize);
								} else {
									$outstr .= "\t(s)\n";
								}
							} else {  ## only mark for writing/saving (global Max Number set)
								if(@doneIntervs < $gFl) {
									$outstr .= "\t(marked for saving)\n";

									push(@{$allRegionsBySc{$QNkey}->{$scKey}},  [ ($regNaame
													, $entr
													, $cchr
													, $distribution->{$entr->[0]}->[5]
													, @doneIntervs +1
													) ] );
								} else {
									$outstr .= "\n";
								}
							}
						} else {
							$outstr .= "\n";
						}

						push(@doneIntervs, [ ($entr->[1], $entr->[2]) ] );
						if((! $writeFlag) && (@doneIntervs >= $nFl)) {
							++$breakFl;
							last;
						}
					}
					last if($breakFl);
				}
				last if($breakFl);
			  }
			  last if($breakFl);
		        }
			if(@doneIntervs) {
				print "chromosome ${cchr}:\n" . $outstr ;
			} else {
				$|=1;
				print ".";
			}
			## go to next chr
		}
       }
       	print "files summary:\n${savedFiles}" if($savedFiles);

       if($gFl) { ## write out globally max gFl best hits (by score)  ## have to add ordering by QN first!!!
		my $breakFl=0;
		my $doneRegs=0;
		my $outstr="";
		my $savedFiles="";
	       #foreach my $scKey (reverse (sort {$a <=> $b} (keys %allRegionsBySc) ) ) {
	       foreach my $QNkey (reverse (sort {$a <=> $b} (keys %allRegionsBySc) ) ) {
			foreach my $scKey (reverse (sort {$a <=> $b} (keys %{$allRegionsBySc{$QNkey}}) ) ) {
			 foreach my $entr  (@{$allRegionsBySc{$QNkey}->{$scKey}}) {
				my $regNaame=$entr->[2] . "\t" . $entr->[0];
				$outstr .= $regNaame;
				if($writeFlag && ($doneRegs < $gFl)) {  ### save region to syntenyRegion-File...
					$outstr .= "\t(saved)\n";
					# write region to file
					$savedFiles .= writeSynRegToFile($entr->[4],  $regNaame, $entr->[3], $entr->[1], $entr->[2],  $realSlength, $flP, $doneRegs +1, $flSize);
				} else {
					$outstr .= "\n";
				}
				++$doneRegs;
				if((! $writeFlag) && ($doneRegs >= $gFl)) {
					++$breakFl;
					last;
				}
			}
			last if($breakFl);
		}
			last if($breakFl);
	}

		 if($doneRegs) {
			if($restrChr) {
				print "\nchr_${restrChr}-best-${gFl}-regions:\n=============================\n";
			} else {
				print "\nglobally-best-${gFl}-regions:\n=============================\n";
			}
			print "chr\tSumOfMaxScores\tstart\trstart\tend\tQN\ttotalHSPcount\ttotalscore\t%ofRange\n" . $outstr;
		} else {
			print "Nothing detected!\n";
		}
		print "global files summary:\n${savedFiles}" if($savedFiles);
	}
} else {  # scOrderFlag set (i.e. old version without taking into account any score information per se! Not recommended actually anymore!
	#
	# obsolete version!!!
	#
	#
	print "chromosome\tregionID\tsaveFlag\n############################\n";
	foreach  my $cchr (keys %blastHits) {
		next if($restrChr && ($restrChr ne $cchr));

		my $rToSave=$nFl +1; # count down regions to save for current chr

    # approach is: find syntenyRegions and start with restriction to find all queries
    # down to only the specified percentage of queries within regions
    # should prevent double information and overlapping regions
    my $cQN = $exonNumb +1; # start to find all-query-regions...
    my @sortedRegionsL = (); # to be a sorted list (sort index queryNs, followed by sort by length) of keys to blasthits in the blastHits-hash
    my %countRegionsH = ();
   ##########
    # get to know upper limit where to start with search
    my @normalRHversions=();
    my $res=1;
    my $preTqn=$minQN;
    while($res && ($preTqn < $cQN)) {
		 my %regionsHash=();
  	 filterHash2SyntenicRKeys2(\%regionsHash, $blastHits{$cchr}, \@exonList, $flSize, ${preTqn});
		 if($res=getHashTotalEntries(\%regionsHash, 0)) {
	     push(@normalRHversions, \%regionsHash);
	   }
	   ++${preTqn};
	   print STDERR ".";
    }
    $cQN = ${minQN} + int(@normalRHversions); ### now maximum available QN determined
    ################
    #    print "\nstarting procedure for chromosome ${cchr}...\n";
        ###print "contig/chromosome ${cchr}...\n";
    my $nCounter=@normalRHversions;
    #print "cQN=${cQN}, minQN=${minQN}, nCounter=${nCounter}\n";
    while(${cQN} > ${minQN}) {
	  #print "cQN=${cQN}, minQN=${minQN}, nCounter=${nCounter}\n";
	   --${cQN};
	   --$nCounter;
	  	if($countRegionsH{$cQN}) {
	     $normalRHversions[$nCounter]=undef;
	     next;
	    }
	  	#print "START with cQN=${cQN}...\n";
    	my %regionsHash= %{$normalRHversions[$nCounter]};
		### new idea: rank the regions acc to diffQNumber and only recalculate for parts that
		### overlap within Qnumber-group, if wished.
      if($res=getHashTotalEntries(\%regionsHash, 0)) {
	    # order (non-overlapping) results according to numberOfDifferentQueryHits (the more the better)
	    # and synRegion-length (the shorter the better, when equal QueryHits)
	      my %sortedRKeys = ();
	      groupRegions(\%sortedRKeys, \%regionsHash);
	      for my $synQN (reverse sort {$a <=> $b} (keys %sortedRKeys)) {
		      next if( exists($countRegionsH{$synQN}));
		      my @cRList=();
		      for my $synL (sort {$a <=> $b} (keys %{$sortedRKeys{$synQN}})) {
		       foreach my $ck (@{$sortedRKeys{$synQN}->{$synL}}) {
			      push(@cRList, $ck);
		       }
		      }
		      if(! detectRegionOverlap(\@cRList)) {
		      # current group (acc to diffQnumber) does not overlap
		       foreach my $ck (@cRList) {
			      if(! detectRegionOverlap( [ ($ck, @sortedRegionsL) ] )) {
			       # current group does not overlap with previous saved results
			       push(@sortedRegionsL, $ck);
			       ### print $ck;
			       ++$countRegionsH{$synQN};
			       if($writeFlag && ((--${rToSave}) > 0)) {
				       # write region to file
				       writeSynRegionToFile($ck, \%regionsHash, $cchr, $nFl - $rToSave +1);
				       ###print " (saved)";
			       }
			       ###print "\n";
			      }	else {
			         print "region \"$ck\" does overlap with a previous result -> skipped\n";
			      }
		       }
		      } else {
		        if($cQN == $minQN) { # smallest allowed queryNumber-loop reached?
			    		print "Some entries on diffQ-Level \"${synQN}\" do overlap.\n";
			        my $startQN=${synQN};
			        #print "skipped at this moment...\n";
			        my $iflP=$flP;
			        while(detectRegionOverlap(\@cRList) && ($iflP > 0.1) ) {
			         print "\$iflP decreased to " . ($iflP=(int($iflP*10)-1)/10.0) . "\t";
			         %regionsHash=();
			         #filterHash2SyntenicRKeys2(\%regionsHash, $blastHits{$cchr}, \@exonList, $iflP * $slength, $cQN);
			         filterHash2SyntenicRKeys2(\%regionsHash, $blastHits{$cchr}, \@exonList, $iflP * $realSlength, $cQN);
			         if($res=getHashTotalEntries(\%regionsHash, 0)) {
				       	my %newSortedRKs=();
				        groupRegions(\%newSortedRKs, \%regionsHash);
				        foreach my $sQN (reverse sort {$a <=> $b} (keys %newSortedRKs)) {
				         next if(($sQN > $startQN) || exists($countRegionsH{$synQN}));
				         @cRList=();
				         for my $sL (sort {$a <=> $b} (keys %{$newSortedRKs{$sQN}})) {
					        foreach my $ck (@{$newSortedRKs{$sQN}->{$sL}}) {
					         push(@cRList, $ck);
					        }
				         }
				         if(! detectRegionOverlap(\@cRList)) {
					        # current group (acc to diffQnumber) does not overlap
					        foreach my $ck (@cRList) {
					         if(! detectRegionOverlap( [ ($ck, @sortedRegionsL) ] )) {
						     # current group does not overlap with previous saved results
						        push(@sortedRegionsL, $ck);
						        ###print $ck;
						        ++$countRegionsH{$sQN};
						        if($writeFlag && ((--${rToSave}) > 0)) {
						        # write region to file
						         writeSynRegionToFile($ck, \%regionsHash, $cchr, $nFl - $rToSave +1);
						         ###print " (saved)";
						        }
						          ###print "\n";
					         } else {
						          print "region \"$ck\" does overlap with a previous result -> skipped\n";
					         }
					        }
				         } else {
					         print "Some entries on diffQ-Level \"${sQN}\" do still overlap.\n";
					         $startQN=${sQN};
					         last; # do next iteration with decreasing iflP... while-loop
				         }
			        } #foreach my $sQN
			    }  #if res...
			   } ## while detectRegionOverlap downwards...

			   if(detectRegionOverlap(\@cRList)) {
			    print "Error: Could not eliminate overlapping by decreasing SizeFactor !\n"
				      , "Aborted at diffQ-Level \"${startQN}\"\n";
			   }
		    } #if($cQN == $minQN)
		    last;
		   } # else  of (if(! detectRegionOverlap(\@cRList)))
	    } #for my $synQN
	   } # if($res=getHashTotalEntries(\%regionsHash...
   	$normalRHversions[$nCounter]=undef;
   }
   if(@sortedRegionsL) {
	 	for(my $synIdx=0; $synIdx<@sortedRegionsL; ++$synIdx) {
	    print "${cchr}\t" . $sortedRegionsL[$synIdx] . (($writeFlag && ($synIdx < $nFl)) ? ("\t(saved)\n") : "\n");
	  }
   }
} #foreach  my $cchr (keys %blastHits...

}


print "\n# Done!\n";
exit 0;

#################################################################################################

sub HELP_MESSAGE {
	pod2usage(-verbose => 1);
}
sub VERSION_MESSAGE {
	print $main::VERSION;
}

################################################################################

sub calcSlidingWindowStat {
	my($blastHits, $flP, $refSize, $minQN, $mScThreshold, $breadth, $stepSfrac)=@_;
	### blastHits are in hash-ref format of  filterBlastFile2HashEq , scoreFl defines, if the bitscore is used instead of simple counts

	$mScThreshold=0 unless(defined($mScThreshold));
	#my $breadth= $flP * $refSize;

	$stepSfrac = 100. unless(defined($stepSfrac));

	my @sortedHashKs=sort keys %{$blastHits};



	my $globalStart=SynBlast::MyUtils::max(int(getKeyPart($sortedHashKs[0], 1))-$breadth + 1, 1);
	my $globalEnd=SynBlast::MyUtils::max(int(getKeyPart($sortedHashKs[-1], 1))-$breadth + 1, 1);
	my %distribution=();


	#my %dummyC=();
	#for(my $cHitIdx=0; $cHitIdx<@sortedHashKs; ++$cHitIdx) {
	#	++$dummyC{getKeyPart($sortedHashKs[$cHitIdx], 5)};
	#}
	#return(undef) unless(scalar(keys %dummyC) >= $minQN);

	return(undef) unless($globalStart <= $globalEnd);

	my $contSkips=0; ## to skip some failing entries later on
	################################################# step-size critical
	my $stepSize=int($refSize/$stepSfrac)+1;

	##
	#print STDERR "\n\nGSTART= ${globalStart}, GEND= ${globalEnd}, RefSize=${refSize}, StepSize= ${stepSize}, $breadth= ${breadth}\n\n";
	##

	for(my $iii=$globalStart; $iii<=$globalEnd   + $stepSize; $iii+=$stepSize) {
	     my %intervalC=();  ## number of hits per query
	     #my %intervalSc=(); ## sum of scores per query
	     my %intervalMSc=(); ## keeps max-score per query
	     my $curSkips=0;
	     my $matchFl=0;
	     my $lastChitIdx=$contSkips;
	     my $realStartPos=undef;  ## contains real start pos (smallest start pos used within region $iii)

	     for(my $cHitIdx=$contSkips; $cHitIdx<@sortedHashKs; ++$cHitIdx) {
		     my $currHitStartP=getKeyPart($sortedHashKs[$cHitIdx], 2);
			if( ( $currHitStartP >= $iii) && (getKeyPart($sortedHashKs[$cHitIdx], 1) <= $iii + $breadth)) { ### curBlastHit fully contained in interval?
				my $qry= getKeyPart($sortedHashKs[$cHitIdx], 5);
				my $scoore=int($blastHits->{$sortedHashKs[$cHitIdx]}->[-1]);


					#if(! exists($intervalC{$qry})) {  #no hit of same queryID already counted for current interval?
					++$distribution{$iii}->[0];  # count total sum of HSPs (hits of any queryID)
					#}
					# sum up bit-score values of hits of same interval (only considered in sum if score above a threshold)
					if($scoore >= $mScThreshold) {
						++$intervalC{$qry};
						$distribution{$iii}->[1] += $scoore;
						#$intervalSc{$qry}+= $scoore;
						$intervalMSc{$qry} = SynBlast::MyUtils::max($scoore, $intervalMSc{$qry});
					}

				++$matchFl;
				push(@{$distribution{$iii}->[5]}, $blastHits->{$sortedHashKs[$cHitIdx]});

				$lastChitIdx=$cHitIdx;
				$realStartPos = SynBlast::MyUtils::min($currHitStartP, $realStartPos);

			} elsif(!  (getKeyPart($sortedHashKs[$cHitIdx], 1) < $iii + $breadth) ) { ### curBlastHit not fully contained in interval due to intervalEnd condition
						###-> all following hits will fail too
						### so skip loop already here
				last;
			} elsif(! (getKeyPart($sortedHashKs[$cHitIdx], 2) >= $iii)) {
				### count number of continous fails from beginning (minus previously accounted skips)
				###not clear####
				++$curSkips if(! $matchFl);
			}

		}
	    ## here now all values for current interval start calculated...
	    #print STDERR "slidingWindow[" . join(",", ($iii, $contSkips, $lastChitIdx, $globalStart, $globalEnd))  . "]\n";
	    if(! defined($distribution{$iii}->[1])) {
		    $distribution{$iii}->[0]=0;
		    $distribution{$iii}->[1]=0;
		    $distribution{$iii}->[2]=0;
		    $distribution{$iii}->[3]=0;
		    $distribution{$iii}->[4]=0;
		    $distribution{$iii}->[5]= [()];
		    $distribution{$iii}->[6]=0;

	    } else {
		    foreach my $quuery (keys %intervalMSc) {
			    #$distribution{$iii}->[4] += (int(10 * $intervalSc{$quuery} / $intervalC{$quuery})/10.);   ## add score-per-query value to a "weighted sum of scores"
			    $distribution{$iii}->[4] += $intervalMSc{$quuery}  # now use: sum of max-scores-per-query in region
		    }
	    $distribution{$iii}->[2]= scalar(keys %intervalC); # save also the number of different queryIDs... (that have at least one representative above threshold
	    $distribution{$iii}->[3]= int(getKeyPart($sortedHashKs[$lastChitIdx], 1)); # save also the genomic end position really used of/for current interval
	    $distribution{$iii}->[6]= $realStartPos;  # save also the genomic start position really used of/for current interval
		}
	     $contSkips += int($curSkips/2.);  ## add skips from current/last session to global skipCount
						## only half due to inexact/larger stepSizes
    	}


	return(\%distribution);

}


sub writeSynRegToFile {
	my($regNr, $regNaame, $regHits, $entry, $mychr, $realS, $flP, $globRegNr, $flSize)=@_;
	## entry is of format
	#start\trealstart\tend\tQN\ttotalHSPcount\tscore\t%ofRange
	$regNr = sprintf("%02d", $regNr);
	my $filestr = "${DBname}.${mychr}.minQN${minQN}.fac${flP}.${regNr}.blastresult";
	if(defined($globRegNr)) {
		$globRegNr = sprintf("%02d", $globRegNr);
		#$filestr =  "${DBname}.glob${globRegNr}.${mychr}.syn${minQP}.fac${flP}.${regNr}.blastresult";
		$filestr =  "${DBname}.glob${globRegNr}.${mychr}.minQN${minQN}.fac${flP}.${regNr}.blastresult";
	}
	##print "writing regions to file
	open(FBLASTF, '>', "${synResultDir}/${filestr}")  or die "Error while creating \"${synResultDir}/${filestr}\": $!\n";
	print FBLASTF "# ${regNaame}\n"
				  . "# ${filestr}\n"
				  . "# generated by: " . $cmdlineWas . "\n"
				  ."##\n";
	# print FBLASTF "#ensemblDate " . $EnsemblD
	# 				. "\n#synOrg " . $DBname
	# 				. "\n#synChr " . $mychr . "\n";
  	print FBLASTF "#synStart " . $entry->[1]
  				. "\n#synEnd " . $entry->[2]
  				. "\n#queryCount " . $entry->[3]
				. "\n#HSPcount " . $entry->[4]
  				. "\n#maxExtension ${flP}x${realS} (" . ($flSize)  . ")"
				. "\n#PercentUsed " . $entry->[6] . "\n##\n";
	foreach my $blhiit (@{$regHits}) {
		print FBLASTF join("\t", @$blhiit) . "\n";
	}
	close FBLASTF;
	return($filestr . "\n");
}

################################################################################
# reads all blast-data into hash of hashes (grouped by chromosome and position_data)
# optionally filtered on certain value at certain idx-position
# wublast-flag to read also wublast-tabular format
sub filterBlastFile2HashEq {
    my ($hashref, $filestr, $idx, $value, $wuformatFl) = @_;

    my $queryIdx = 6;
    my $targetIdx = 8;

    if($wuformatFl) {
	    $queryIdx = 6; #8;
	    $targetIdx = 8; #10;
    }

    my $queryIdxs={()};
    $queryIdxs = $idx if(defined($idx) && (! defined($value)));  ### if only 3rd param set, this is the hash of allowed querys

    open(BLASTFILE, '<', "${filestr}")  or die "Error while reading \"${filestr}\": $!\n";
    while(<BLASTFILE>) {
	if ($_ =~ m/^[^\#]/) {
	    chomp;
	    my @curHit=split; #(/\s+|\t/, $_);
	    #process curHit to get matches in chromosome orientation

	    	  if($wuformatFl) {  ### 22 strd fields, optionally two more (24)
		  next unless((@curHit >= 22) && (@curHit <= 24));
		  # Fields: qid   sid     E       N       Sprime  S       alignlen        nident  npos    nmism   pcident pcpos   qgaps   qgaplen sgaps   sgaplen qframe  qstart  qend    sframe  sstart  send    group   links
		   @curHit = ( $curHit[0]  # pos0: queryID
					,$curHit[1] #pos1: subjectID
					,$curHit[10]  #pos2:perc.id
				#	,$curHit[11] #pos3: perc.pos
				#	,$curHit[16] . "/" . $curHit[19] #pos4: qry/sbj-frame
					,$curHit[6]  #pos5: aln-length
					,$curHit[9]   #pos6: nmism
					,$curHit[14]  #pos7: gap opens?
					,$curHit[17]  #pos8:q.start
					,$curHit[18]   #pos9:q.end
					,$curHit[20]   #pos10:s.start
					,$curHit[21]   #pos11:s.end
					,$curHit[2]    #pos12:evalue
					,$curHit[4]     #pos13:bitscore
				);
		}



	    next unless( (! scalar(keys %{$queryIdxs})) || exists($queryIdxs->{$curHit[0]}));   ## skip this entry if queryIdxs is defined and does not contain current queryID

	    do { SynBlast::MyUtils::swapArray(\@curHit, $queryIdx, $queryIdx+1);
			SynBlast::MyUtils::swapArray(\@curHit, $targetIdx, $targetIdx+1);
		} unless ( int($curHit[$targetIdx+1]) >= int($curHit[$targetIdx]));




	    my $curH9=SynBlast::MyUtils::getLeadingZeroStr($curHit[$targetIdx+1],10);
#qx(printf "%.10d" $curHit[9]);
	    # my $curH8=qx(printf "%.10d" $curHit[8]);
	    # my $curKeyChrHits="${curH9}_${curHit[8]}_${curHit[0]}";

	    my $curKeyChrHits="${curH9}_${curHit[$targetIdx]}_${curHit[$queryIdx]}_${curHit[$queryIdx +1]}_${curHit[0]}";
	    ###my $curKeyChrHits="${curH9}_${curHit[8]}_0_0_${curHit[0]}";


	    #if ( (defined($value) && ("$curHit[$idx]" eq "${value}")) || (! defined($value)) ) {
	    if ( (! defined($value)) || ("$curHit[$idx]" eq "${value}") ) {
		if (defined($hashref->{$curHit[1]}->{${curKeyChrHits}})) {
		    #print STDERR "ERROR???!!!! Updating:\n";
		    #print STDERR "OldEntry= " . join(" ", @{ $hashref->{$curHit[1]}->{${curKeyChrHits}} } ) . "\n";
		    #print STDERR "NewEntry= " . join(" ", @curHit) . "\n";
		  print STDERR "E";
		}
		$hashref->{$curHit[1]}->{${curKeyChrHits}} = [ @curHit ] ;
	    }
	}
    }
    close(BLASTFILE);
}


sub getKeyPart {
    my ($keystr, $pos) = @_;
    if ($keystr =~ m/^(\d+)_(\d+)_(\d+)_(\d+)_(.*)$/) {
	return($1) if ($pos == 1);
	return($2) if ($pos == 2);
	return($3) if ($pos == 3);
	return($4) if ($pos == 4);
	return($5) if ($pos == 5);
    }
}





#####################################################################################
#### old subs (for obsolete filter version)

sub groupRegions {
    my ($outHRef, $regionH) = @_;
    %{$outHRef} = ();
    foreach my $synR (keys %{$regionH}) {
	my @synR = split(/_/, $synR);
	my $cQueryN = $synR[2];
	my $cRlength = abs($synR[1] - $synR[0]) + 1;
	push(@{$outHRef->{$cQueryN}->{$cRlength}}, $synR);
    }
}


sub writeSynRegionToFile {
    my ($synRidx, $regHRef, $cChrom, $cFileN) = @_;
    #my $filestr = "${synResultDir}/${cChrom}.${synRidx}.synteny_${minQP}.${DBname}.blastresult";

    my $filestr = "${synResultDir}/${DBname}.${cChrom}.syn${minQP}.${cFileN}.blastresult";
    ##print "writing regions to file
    open(FBLASTF, '>', "${filestr}")  or die "Error while creating \"${filestr}\": $!\n";
    print FBLASTF "# ${synRidx} \n"
    					##	. "# ${cChrom}.${synRidx}.synteny_${minQP}.${DBname}.blastresult\n##\n";
    					#. "# ${DBname}.${cChrom}.syn${minQP}.${cFileN}.blastresult\n##\n";
					. "# ${DBname}.${cChrom}.minQN${minQN}.${cFileN}.blastresult\n##\n";
    print FBLASTF "#ensemblDate " . $EnsemblD
    						. "\n#synOrg " . $DBname
    						. "\n#synChr " . $cChrom . "\n";
  #  print FBLASTF "#synStart " . getKeyPart($synRidx, 1) . "\n";
  #  print FBLASTF "#synEnd " . getKeyPart($synRidx, 2) . "\n";
  	my @csynIdx=split(/_/,   $synRidx);
  	print FBLASTF "#synStart " . $csynIdx[0]
  				. "\n#synEnd " . $csynIdx[1]
  				. "\n#queryCount " . $csynIdx[2]
  				. "\n#maxAllowedSize " . $csynIdx[3] . "\n##\n";

    foreach my $chrRk (@{$regHRef->{$synRidx}}) {
	print FBLASTF join("\t", @{ $blastHits{$cChrom}->{$chrRk} } ) . "\n";
    }
    close FBLASTF;
    print "\"${filestr}\" created.\n\n";
}


sub detectRegionOverlap {
    my ($regKeysRef) = @_;
    return(0) if(! @{$regKeysRef});
    my $ovFound = 0;
    for(my $lii = 0; ($lii < (@{$regKeysRef} -1)) && (! $ovFound); ++$lii) {
	my @iRkey= split(/_/, $regKeysRef->[$lii]);
	my $iRstart= $iRkey[0];
	my $iRend= $iRkey[1];
	for(my $ljj = $lii + 1; ($ljj < @{$regKeysRef}) && (! $ovFound); ++$ljj) {
	    my @jRkey= split(/_/, $regKeysRef->[$ljj]);
	    my $jRstart= $jRkey[0];
	    my $jRend= $jRkey[1];
	    $ovFound += (SynBlast::MyUtils::getOverlapSize($iRstart, $iRend, $jRstart, $jRend) > 0);
	}
    }

    return($ovFound);

}


sub checkKeyIsInside {
    my ($keystr, $lborder, $rborder) = @_;
    if ($keystr =~ m/^(\d+)_(\d+)_(.*)$/) {
	return((int($1) <= $rborder) && (int($2) >= $lborder));
    }
}


################# old current version ##################
################################################
sub filterHash2SyntenicRKeys2 {
    my ($href, $hashref, $querylist, $maxSize, $minQN) = @_;

    if (! defined($minQN)) {
	   $minQN=2;
    }
    my $querystr=join(" ", @$querylist);
    return(undef) if (@$querylist < 1);

    my @sortedHashKs=sort keys %{$hashref};
    my $curStartIdx=@sortedHashKs -1;
    $maxSize=int($maxSize);
    my $oriMaxS = $maxSize;
    my $lastdbEnd=undef;
    my $lastdbStart=undef;
    #while($curStartIdx >= @{$querylist} - 1) {
	  # as long there are at least minQN entries possible...
	  while($curStartIdx >= $minQN - 1) {
	   my @curKeys=();
	   my @curKeyIdxs=();
	   my %curQueryCount=();
	   while(($curStartIdx >= 0) && (! (index($querystr, getKeyPart($sortedHashKs[$curStartIdx],5)) >= 0)) ) { # no match in queryList?
	    --$curStartIdx;
	   }
	   my $curRend=int(getKeyPart($sortedHashKs[$curStartIdx], 1));
	   my $curRstart=SynBlast::MyUtils::max(1, $curRend - ($maxSize));
	   ### determine smallest index fulfilling regionsize ( blastHit-endpos >= curRstart)
	   my $curEndIdx=evalWindow(\@sortedHashKs, $curStartIdx, $curRstart, $curRend, int($curStartIdx /4)+1);
     my $curRdbmax=$curRend;
	   my $curRdbmin=$curRend;
# go through current index intervall and store protein occurencies
	   for(my $ii=$curEndIdx; $ii <= $curStartIdx; ++$ii) {
	    if( checkKeyIsInside($sortedHashKs[$ii], $curRstart, $curRend)) { # lies in curRegion...
		   my $curQuery=getKeyPart($sortedHashKs[$ii],5);
		   if(index($querystr, $curQuery) >= 0) { # match in queryList?
		    $curQueryCount{$curQuery}++;
		    push(@curKeyIdxs, $ii);
		    push(@curKeys, $sortedHashKs[$ii]);
		    $curRdbmin=SynBlast::MyUtils::min($curRdbmin, int(getKeyPart($sortedHashKs[$ii],2)));
		   }
	    }
	   }
	   my $differentQs=keys(%curQueryCount);
	my $rkey="${curRdbmin}_${curRdbmax}_${differentQs}_${maxSize}_${minQN}";
	#print "number of entries in current window: " . @curKeyIdxs . " bzw. " . @curKeys . "\n";
	if( ${differentQs} >= $minQN) {
	    # number of different queries is at least $minQP percent of total query number?
	    # save curRegion-data and keys in output-referencehash...
	    if (defined($href->{$rkey})) {
		print "Error?: Same Region with different starting point reached: ${rkey}\n";
		print "Number of Entries was: " . @{$href->{$rkey}} . "\n";
		print "New Number of Entries is: " . @curKeys . "\n";
	    }
	    if(($lastdbEnd) && (($lastdbStart == $curRdbmin) && ($curRdbmax < $lastdbEnd))) {
		#only a subregion of last one found -> do not save ...
		#print "obmitting subregion...\n";
	    }
	    else {
		$href->{$rkey} = \@curKeys;
		$lastdbEnd=${curRdbmax};
		$lastdbStart=${curRdbmin};
	    }
	}
	#else {
	    # number of different queries is lower than $minQP percent of total query number?
	#}
	if(! defined($curKeyIdxs[0])) {
	    $curStartIdx--;
	}
	else {
	    my $smHitIdx=$curKeyIdxs[0];
	    my $priorSHitEnd= int(getKeyPart($sortedHashKs[$smHitIdx -1], 1));
	    if($smHitIdx > 0) {
	    	if( ${differentQs} == 1) { # only one match type -> go to last match
		     if($curStartIdx==$smHitIdx) {
			    $curStartIdx--;
		     } else { $curStartIdx=$smHitIdx; }
		    } else { # more than one matching query name
		      if(($curRdbmin - $priorSHitEnd) > int($maxSize)) {
			# if space between region and next putative blastHit is to large(?), proceed directly from the next putative region-start...
			     $curStartIdx=$smHitIdx -1;
			     $maxSize=$oriMaxS;
		      }
		      else {
			     while(($curStartIdx >= 0) && $curRend - int(getKeyPart($sortedHashKs[$curStartIdx], 1)) < 1) { #$flankingSize ) {
			      $curStartIdx--;
			     }
			    }
		    }
	    } else {
		   $curStartIdx=-1;
	   }
	}
 } #while($curStartIdx >= @{$querylist} - 1)...
}

sub evalWindow { # returns lowest Idx of Key-Array-Entry whose blastHitEnd is still at least startpos and hence a candidate for curWindow-Hits
    my ($keysref, $startidx, $startpos, $endpos, $stepsize) = @_;
    $stepsize = 1 if(! defined($stepsize));
    my $ii=$startidx;
    while(($ii >= 0) && (getKeyPart($keysref->[$ii], 1) >= $startpos)) { # while blastHitEnd is larger than startpos
	$ii = $ii - $stepsize;
    }
    return(evalWindow($keysref, $ii + $stepsize, $startpos, $endpos, int($stepsize/2))) if ($stepsize > 1);
    return($ii + $stepsize);
}


sub getHashTotalEntries {
    my ($hashref, $printFlag) = @_;
    my $curE=keys(%$hashref);
    if ($printFlag) { print $curE . "\n"; }
    return($curE);
}

sub getListTotalEntries {
    my ($listref, $printFlag) = @_;
    my $curE=@$listref;
    if ($printFlag) { print $curE . "\n"; }
    return($curE);
}
#### end of obsolete subs
############################################################################################################
