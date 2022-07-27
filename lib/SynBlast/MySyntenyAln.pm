#!/usr/bin/env perl
# Last changed Time-stamp: <2008-08-19 15:05:22 joe>
# 2006-03-04
# joe@bioinf.uni-leipzig.de
########################################################################

use strict;
use warnings;

package SynBlast::MySyntenyAln; 
require Exporter;
my @ISA = qw(Exporter);
my @EXPORT = qw(getSyntenyAln getBestBlastLoci);
my @EXPORT_OK = @EXPORT;

use SynBlast::MyUtils;
use SynBlast::MySyntenyUtils;
use SynBlast::EnsEMBLaccess;

use GD;

##########################################################################


##########################################################################
# chain blast hit loci of blastFile, main part used by getBestBlastHits.pl
##########################################################################
sub  getBestBlastLoci {
  my (  $infoDataR    ## list of reference region infos
	, $blastFile 	## blast hits in tabular format on local disk
	, $blastQID 	## queryID to determine the best loci for
	, $scutoff   	## score cutoff (chained loci must reach at least this threshold of bitscore sum)
	, $protOvl		## maximal allowed aa-overlap size between assigned HSP intervals in the chain
	, $verboseFl
	, $restrChr        ## restrict to a certain chromosome/contig?
	, $nbest		## only report the n-best regions
	, $comparaFl    ## retrieve compara information?
	### , $comparaInitF
	, $organismN    ## target species/db name
	, $QgeneID	
	, $QgeneName
	, $rangeFact	## gene locus range factor
	, $deepFl		## use old version for chaining procedure?
	, $sliceAdapt    
	, $comparaData  ## pre-stored comparaData
	, $globMaxSize   ## set a maximal size globally for each locus
	, $wublastFl ## expects then wublast tabular format instead of ncbi format
	)= @_;
  
    
    my ($queryStI, $targetStI)=(6, 8);
    ($queryStI, $targetStI)=(8, 10) if($wublastFl);
       
    ### to use old version of calcSplAln (deep recursion, much slower)
    $deepFl=0 if(! defined($deepFl));
    
    $globMaxSize=0 if(! defined($globMaxSize));
    
    $rangeFact=2 if(! defined($rangeFact));
    $verboseFl=1 if(! defined($verboseFl));
    $comparaFl=0 if(! defined($comparaFl));
    #$sFl=0 if(! defined($sFl));

    $nbest = 1 if(! defined($nbest));
    
    $protOvl=20 if(! defined($protOvl));
    my %protHits=();

    
    # here use chromosome/contig name as grouping index (1) instead of query-name (0)
    if(! (SynBlast::MySyntenyUtils::blastfile2sortedLists("${blastFile}"
													, \%protHits
													, 1
													, $targetStI + 1  ## sort idx!!!
													, 0
													, $blastQID
													, undef
													, undef
													, undef
													, undef
													, $wublastFl) > 0)) {
				print STDERR ">warning: ${blastFile} contains no data (at least for query \"${blastQID}\"! (getBestBlastHit)\n";
				return(undef, undef);
	}

     
    # extract infoFileData
   my ($exons, $queryIdx, $refExons, $exonList, $infoH, $icolIdx, $realFlankS, $realSlength) = @{$infoDataR};
     
    ## calculate for each chromosome/contig with hits corresponding loci and overall score of consistent blast-hit series
    ## for proteinQueryID $blastQID
    
    if(! exists($queryIdx->{$blastQID})) {
	warn "WARNING: ${blastQID} not a key in queryIdx (maybe non-compatible infofile and blastresults used?!)\n";
    }	    
    my $refEntry=$exons->[$queryIdx->{$blastQID}];
    my $refLength=abs($refEntry->[$icolIdx->{'endPos'}] - $refEntry->[$icolIdx->{'startPos'}]) +1;
    #print "refLength=", $refLength, "\n";
    
    my $maxsize=$rangeFact * $refLength;
    $maxsize=SynBlast::MyUtils::min($globMaxSize, $maxsize) if($globMaxSize);
    
    my %bestBlastLoci=();

    foreach my $curChrom (sort keys %protHits) {
	next if(defined($restrChr) && ($restrChr ne $curChrom));

	#print "CURRENT CONTIG/CHROMOSOME: ${curChrom}\n-------------------------------------\n";
	my %lociScores=(); # contains later the hash of list of hash of SplAln-scores for pair (queryProt,targetLoci)
	if(! $deepFl) {
	  calcSplAln2(\@{$protHits{$curChrom}}
				, \%lociScores
				, $blastQID
				, $maxsize 
				, $scutoff
				, $protOvl
				, 0
				, $refEntry->[$icolIdx->{'protLength'}]
				, $queryStI
				, $targetStI				
				);#2);#0);
	} else {
	  print STDERR ">using old version for blast-hit assembly...\n";
	  calcSplAln(\@{$protHits{$curChrom}}, \%lociScores, $blastQID, $maxsize , $scutoff, $protOvl, 0, $refEntry->[$icolIdx->{'protLength'}], $queryStI, $targetStI);#2);#0);
	}
	
	my @lkey = keys %lociScores;
	#print STDERR "No hits matching specified constraines could be found! (queryID was ${blastQID} on chr ${curChrom})\n" unless(@lkey);
	next if(! @lkey);
	die(("Error: no more than one protQueryID allowed for blastfile in this appl.! (content of lkey was " . join(";", @lkey) . ")\n")) unless (@lkey <= 1);
	foreach my $cEntryH (@{$lociScores{$lkey[0]}}) {
	    push(@{$bestBlastLoci{int(($cEntryH->{'score'}) * 10)}}, $cEntryH);

	}
    }

    ### sort desc to key 'score'
    my @sortdKeys = reverse sort {$a <=> $b} (keys %bestBlastLoci);
    
    my $noLoci= SynBlast::MyUtils::min($nbest, int(@sortdKeys));
    #my $noLoci=$nbest;
    
    if($verboseFl) {
    	print "#chr\tstartPos\tendPos\tori\thitsCount\thitsLength\tqueryCov\tmeanIdentity\tscore";
      print "\tcompara" if($comparaFl);
      print "\n";
    }
    
    
  my @globalHitList=();  
     
 for(my $li=0; $li< $noLoci; ++$li) {
	foreach my $cE (@{$bestBlastLoci{$sortdKeys[$li]}}) {   
	    my %currE= %{$cE}; ## take all entry descriptors 
	    $currE{'hits'}=undef;
	    $currE{'rank'}= $li + 1;
	     
	    if($verboseFl) {
	      print
		  ${QgeneName} . "/" . ${QgeneID}
		  . "\t"
		 .  $cE->{'chr'} 
	      . "\t"
		  . $cE->{'startPos'}
	      . "\t"
		  . $cE->{'endPos'}
	      . "\t"
		  . $cE->{'ori'}
	      . "\t"
		  . $cE->{'hitsCount'}
	     . "\t"
		 . $cE->{'hitsLength'}
	     . "\t"
		 . $cE->{'queryCov'}
	     . "\t"
		 . $cE->{'meanIdentity'}
	     . "\t"
		 . $cE->{'score'}; #. "\n";
	   }
	    
	    if($comparaFl) {
		
		my $coordType = ""; ###undef;
		##$coordType = "chromosome" if(int($cE->{'chr'}));
		
		(my $comparaR, my $resHash) = SynBlast::EnsEMBLaccess::checkComparaData($comparaData->{${QgeneID}}
																		, $sliceAdapt
																		, ${organismN}
																		, $coordType
																		, $cE->{'chr'}
																		, $cE->{'startPos'}
																		, $cE->{'endPos'}
																		, ${QgeneID}
																		, 0
																		, 1);
 		chomp($comparaR);
		print "\t" . $comparaR if($verboseFl);		
		
		$currE{'compara'}= $resHash;
		
	  }
	  
	    print "\n" if($verboseFl);
	    
	    if($verboseFl >1) {
		print "==================_single_hits_==================\n";
		foreach my $cSH (@{$cE->{'hits'}}) {
		    print join("\t", @{$cSH}) . "\n";
		}
		print "\n";
	    }
	
	push(@globalHitList, \%currE);
	}
  }
	return(\@globalHitList);
}










##########################################################################
# combine reference loci into non-overlapping p-loci, by grouping overlapping transcripts of reference proteins 
# according to start/end or mean HSP positions
##########################################################################
sub getCombPlociRef {
    my ($infoDataR, $verboseFl, $sFl)= @_;
    $verboseFl=1 if(! defined($verboseFl));
    $sFl=0 if(! defined($sFl)); # if set, grouping by meanpos

     my ($exons, $queryIdx, $refExons, $exonList, $infoH, $icolIdx, $realFlankS, $realSlength) = @{$infoDataR};
    
   my $refOrgName=$infoH->{"organism"};
   my $refChrCtg=$infoH->{"seqregion"};
    
   my @slicePos = ($infoH->{"from"}, $infoH->{"to"});

    my @plocis=();
    my %plocis=();
 
 for(my $la=0; $la<@{$exons}; $la++) {
	my %tmpHa=();
	my %exons=();
	foreach my $colu (keys %{$icolIdx}) {
		$exons{$colu} = 	$exons->[$la]->[$icolIdx->{$colu}];
	}
	$exons{'chr'} = $refChrCtg;
	
	if(! $sFl) {
	    %tmpHa = (  ($exons->[$la]->[$icolIdx->{'startPos'}] . "_" . $exons->[$la]->[$icolIdx->{'endPos'}]) =>  \%exons); # same format as for qlocis for use in combineLoci2List
	    push(@{$plocis{$exons{'endPos'}}}, \%tmpHa ); # sort index is end-pos on slice
	} else {
	    %tmpHa = ( ($exons->[$la]->[$icolIdx->{'meanPos'}] . "_" . $exons->[$la]->[$icolIdx->{'meanPos'}]) =>  \%exons);
	    push(@{$plocis{$exons{'meanPos'}}}, \%tmpHa ); # sort index is mean-cds-pos on slice (sFlag set)
	}
    }
    
    foreach my $endPos (sort {$a <=> $b} (keys %plocis)) {
	foreach my $haref (@{ $plocis{$endPos} }) {
	    push(@plocis, { %{$haref} });
	}
    }
    %plocis=(); # not needed any more
    
    my @combPlocis=();
    combineLoci2List(\@plocis, \@combPlocis, $sFl, 0);
    
    @combPlocis= reverse @combPlocis; ## to get again increasing order...
    
    if($verboseFl) {
	print "combined p-loci (queries) as follows:\n################################################\n";
	for(my $la=0; $la<@combPlocis; $la++) {
	    print $la . "\t";
	    SynBlast::MyUtils::printHashHashListHash($combPlocis[$la], "=>", "\n", "\n", "\n\n");
	    print "\n\n";
	} 
	
	print "number of query blocks (p-loci) that could be combined into overlapping ones is " . (@plocis - @combPlocis) 
	    . " (old was " . @plocis . ", new is " . @combPlocis . ")\n";
	
	print "number of reference proteins is " . @{$exons} . ".\n";
	
    }
    
    my %qIdxCombPlocis=();
    my @plociDistances=();
    for(my $qIdx=0; $qIdx <@combPlocis; ++$qIdx) {
	my @plkey = keys (%{$combPlocis[$qIdx]});
	(@plkey == 1) or die "Error - There is more than one Hash-Entry for a combined_P_locus (getSyntenyAln)!\n";
	
	if($qIdx) {
	  my @plkey0 = keys (%{$combPlocis[$qIdx-1]});
	  (@plkey0 == 1) or die "Error - There is more than one Hash-Entry for a combined_P_locus (getSyntenyAln)!\n";
	  my($k1from,$k1to)=split(/_/, $plkey[0]);
	  my($k0from,$k0to)=split(/_/, $plkey0[0]);
	  my $kdist= $k1from - $k0to;
	  push(@plociDistances, $kdist);
	}
	
	foreach my $ploc (keys %{$combPlocis[$qIdx]->{$plkey[0]}}) {
	    foreach my $cPE (@{$combPlocis[$qIdx]->{$plkey[0]}->{$ploc}}) { # cPE now a hashref
		$qIdxCombPlocis{($cPE->{'protID'})} = $qIdx;
	    }
	}
    }
    
    print "\n\nReference-Exons/Proteins are/is:\n################################################\n" if($verboseFl);
    my @refPloci=();
    my $refProtStr="";
    my @refLidxs=();
    for(my $li=0; $li<@{$refExons}; $li++) {
				my $refLociN=($qIdxCombPlocis{$refExons->[$li]} +1);
				my $ploci= "p" . $refLociN;
       	my $dummyStr= $refExons->[$li] . " (contained in " . $ploci . ")";
	      $refProtStr .= $dummyStr;
	      $refProtStr .= "; " if($li< (@{$refExons} -1));
  	  if($verboseFl) {
	     print $dummyStr;
	     print ", " if($li< (@{$refExons} -1));
	    }
	    push(@refPloci, [ ($ploci, $refExons->[$li]) ] );
	    push(@refLidxs, $refLociN)
    }
    
    return(\@combPlocis, \@refPloci, $refProtStr, \@refLidxs, SynBlast::MyUtils::max(@plociDistances));
    
}


##########################################################################
# evaluation of blast-hit regions via gene order alignment
# main sub for doSyntenyAlignment.pl, generating all the data and graphics into memory
##########################################################################
sub getSyntenyAln {
    my ($combPlocis
	, $refPloci
	, $refProtStr
	, $refLidxs
	, $infoDataR
	, $blastFile
	, $scutoff, $protOvl, $locSizeFac
	, $verboseFl, $selfcFl, $sFl
	, $pfields, $qfields
	, $gapPp, $gapPq, $mismV, $penaltyFa
	, $regID
	, $ppfields
	, $dpColors
	, $printcFl
	, $compAllFl
	, $initF
	, $compData
	, $blHitsubDir
	, $gRankData
	, $relScFl
	, $deepFl
	, $maxDist
	, $globLocSize
	, $refscores
	)= @_;
    
    $printcFl=0 if(! defined($printcFl));
    $printcFl=1 if($printcFl);

   $verboseFl=1 if(! defined($verboseFl));
    $sFl=0 if(! defined($sFl));
    $selfcFl = 0 if(! defined($selfcFl)); ## flag to save/update best (self)score into reference file
    $ppfields = $pfields if(! defined($ppfields));
    
    $relScFl= (! $selfcFl) if(! defined($relScFl));
    
    $protOvl=30 if(! defined($protOvl));
    $locSizeFac = 2 if(! defined($locSizeFac));
    
    $globLocSize=0 if (! defined($globLocSize));
    
    ### to use old version of calcSplAln (deep recursion, much slower)
    $deepFl=0 if(! defined($deepFl));
    
    # extract infoFileData
    my ($exons, $queryIdx, $refExons, $exonList, $infoH, $icolIdx, $realFlankS, $realSlength) = @{$infoDataR};
    my $refOrgName=$infoH->{"organism"};
    my $refChrCtg=$infoH->{"seqregion"};
    my @slicePos = ($infoH->{"from"}, $infoH->{"to"});
   
    my $ensemblQstr = undef;
    $ensemblQstr = $infoH->{"ensemblRelease"} if(defined($infoH->{"ensemblRelease"}));
########################################

    my($ensemblStr, $synOrg, $synChr, $synStart, $synEnd) = SynBlast::MySyntenyUtils::getSynRegInfo($blastFile);
    
########################################
    my %protHits=();
    # blast hits for each protQuery now explicitely sorted ascending in dna-end-position value!!! (needed for calcSplAln)
    if(! (SynBlast::MySyntenyUtils::blastfile2sortedLists("${blastFile}", \%protHits, 0, 9) > 0)) {
	print STDERR "ERROR: ${blastFile} contains no data! (getSyntenyAln)\n";
	return(undef, undef);
    }

    ## calculate for each query-protein corresponding blast loci and overall score of consistent blast-hit series...
    print STDERR ">chaining blast-hits together...";
    my %lociScores=(); # contains later the hash of list of hash of SplAln-scores for pair (queryProt,targetLoci)
    foreach my $curP (sort keys %protHits) {
	if(! exists($queryIdx->{$curP})) {
	    warn "WARNING: ${curP} not a key in queryIdx (maybe non-compatible infofile and blastresults used?!)\n";
	    next;
	}	    
	my $refEntry=$exons->[$queryIdx->{$curP}];
	my $refLength=abs($refEntry->[$icolIdx->{'endPos'}] - $refEntry->[$icolIdx->{'startPos'}]) +1;
	
	my $maxsize=$locSizeFac * $refLength;
	$maxsize=SynBlast::MyUtils::min($globLocSize, $maxsize) if($globLocSize);
	
	if(! $deepFl) {
	  calcSplAln2(\@{$protHits{$curP}}, \%lociScores, $curP, $maxsize, $scutoff, $protOvl, 0, $refEntry->[$icolIdx->{'protLength'}]);
	} else {
	  print STDERR ">using old version for blast-hit assembly! ...\n";
	  calcSplAln(\@{$protHits{$curP}}, \%lociScores, $curP, $maxsize, $scutoff, $protOvl, 0, $refEntry->[$icolIdx->{'protLength'}]);
	}
	
	
	
	
    }
    
    
    ## get sorted list of loci (sorted by DNA-end-pos increasing) for whole blastregion
    my @locis=();
    qlociSort2List(\%lociScores, \@locis, $sFl);
    
    
    # get combined q-Loci (overlapping loci of different proteins are merged, and only one list-entry per range and protQuery allowed (qLflag set))
    my @combLocis=();
    combineLoci2List(\@locis, \@combLocis, $sFl, 1);
    @combLocis= reverse @combLocis; ## to get again increasing order...
    if($verboseFl) {
	print "combined q-loci (target blocks) as follows:\n################################################\n";
	for(my $la=0; $la<@combLocis; $la++) {
	    print $la . "\t";
	    SynBlast::MyUtils::printHashHashListHash($combLocis[$la], "=>", "\n", "\n", "\n\n");
	    print "\n\n";
	} 
	print "number of target blocks (q-loci) that could be combined into overlapping ones is " . (@locis - @combLocis) 
	    . " (old was " . @locis . ", new is " . @combLocis . ")\n\n";
    }
    print STDERR "done\n";

#######################


    print STDERR ">calculating gene order alignments...";
    ## do alignment ala Needleman-Wunsch with endfree-gaps...
    my $cOptLociRef=\@combLocis;
    my $penaltyFl=1;
    my $whichAlnInfo=0; #"normal (increasing chromosome positions)";
    print "\nAlignment one:\n################################################\n" if($verboseFl);
    
    my @synOrg= split(/_/, $synOrg);
    my $myRegDescr = substr($synOrg[0], 0, 1) .  substr($synOrg[1], 0, 3) . "_";
    $myRegDescr .= "chr" if((length($synChr) < 2) || int($synChr));
    $myRegDescr .= $synChr;
    
    my @refOrg= split(/_/, $refOrgName);
    my $myRefDescr = substr($refOrg[0], 0, 1) .  substr($refOrg[1], 0, 3) . "_";
    $myRefDescr .= "chr" if((length($refChrCtg) < 2) || int($refChrCtg));
    $myRefDescr .= $refChrCtg;
    
    my ($alnScore,$alnStrRef, $alnDPref, $alnIdxR, $alnDPpref)=alignLocis(
							      $combPlocis
							      , \@combLocis
							      , $gapPp
							      , $gapPq
							      , $mismV
							      , $verboseFl
							      , 1
							      , $penaltyFa
							      , $regID
							      , $dpColors
							      , $printcFl
							      , undef
							      , $relScFl
							      , $refLidxs
							      , 1 ## create DPlot
							      , $maxDist  ## maximal distance between columns (only hightlighting, if exceeded)
							      , undef
							      , undef
							      , $myRegDescr
							      , $myRefDescr
							      , $refscores
							      );
    
 
    my @revCombLocis = reverse @combLocis;
    if ($verboseFl) {
	print "\n\n"; print "\nAlignment three (reversed loci):\n################################################\n";
      }
    my ($aln2Score,$aln2StrRef, $aln2DPref, $aln2IdxR, $aln2DPpref)=alignLocis(
								  $combPlocis
								  , \@revCombLocis
								  , $gapPp
								  , $gapPq
								  , $mismV
								  , $verboseFl
								  , 2
								  , $penaltyFa
								  , $regID
								  , $dpColors
								  , $printcFl
								  , undef
								  , $relScFl
								  , $refLidxs
								  , 1
								  , $maxDist
								  , undef
								  , undef
								  , $myRegDescr
								  , $myRefDescr
								  , $refscores
								   );
    
    if($aln2Score > $alnScore) {
	print "-->", "Second alignment (reversed loci) is better!\n" if($verboseFl);
	($alnStrRef, $alnScore, $alnDPref, $alnIdxR, $alnDPpref) = ($aln2StrRef, $aln2Score, $aln2DPref, $aln2IdxR, $aln2DPpref);
	$cOptLociRef= \@revCombLocis;
	$penaltyFl=2;
	$whichAlnInfo=1;
    }
    
    print STDERR "done\n";

   
    my $optHitListRef=getOptHitList($combPlocis, $cOptLociRef);
    
    if ($verboseFl) {
	print "\n\n";
	print "List of optimal Query-Hits (q-numbers according to the optimal alignment) as follows:\n################################################\n";
	foreach my $lid (@{$optHitListRef}) {
	    print "###\n";
	    SynBlast::MyUtils::printListList($lid, "\t" , "\n");
	}
    }
    
    # extract info for reference-Exon-matches (for overview table and comparison of different alignments in results)
    my @qloci=();
    for(my $pli=0; $pli<@{$refPloci}; $pli++) { # for every referenceProt do
	my $pNumber=substr($refPloci->[$pli]->[0],1);
	my $curPID=$refPloci->[$pli]->[1];
	
	for(my $li=0; $li< @{$alnStrRef}; $li++) { # go through alnStr
	    if( "$alnStrRef->[$li]->[0]" eq "$refPloci->[$pli]->[0]" ) { # refExon-position in alnstr found...
		my $qNumber=0;
		my $qStr=$alnStrRef->[$li]->[1];
		$qNumber=int(substr($qStr,1)) if($qStr ne "---");
		if($qNumber > 0) {
		    my @curLociKey=(keys %{$cOptLociRef->[($qNumber -1)]});
		    (@curLociKey == 1) or die "Wrong number of keys! (getSyntenyAln)\n";
		    my %curData=();
		    if(defined($cOptLociRef->[($qNumber -1)]->{$curLociKey[0]}->{$curPID}->[0]) ) {
			%curData= %{$cOptLociRef->[($qNumber -1)]->{$curLociKey[0]}->{$curPID}->[0]}; # is hash-ref
		    } else {
			print STDERR "WARNING: Hash-Reference for q${qNumber} and PID ${curPID} was not defined!\n";
		    }
		    
		    my $curProtName = $exons->[$queryIdx->{$curPID}]->[$icolIdx->{'protName'}];
		    $curData{'protInfo'} = "(${qStr}=>${curPID}" . ( ($curProtName ne "-" ) ? ("=>" . $curProtName ) : "" ) . ")";
		    
		    my $optListEntries = "#For comparison: OptHitEntry for ";
		    for(my $listi=0; $listi< @{$optHitListRef->[$pNumber -1]}; $listi++) {
			if(index($optHitListRef->[$pNumber -1]->[$listi]->[0], $curPID) >= 0) {
			    $optListEntries = $optListEntries . join(" ", @{$optHitListRef->[$pNumber -1]->[$listi]});
			    $optListEntries = $optListEntries . "\n" if($listi< @{$optHitListRef->[$pNumber -1]} -1);
			}
		    }
		    
		    $curData{'OptHitInfo'} = $optListEntries;
		    
		    $curData{'refGeneID'} = getPlocDat($combPlocis->[$pNumber - 1], 0, "geneID");
		    
		    push(@qloci, \%curData); 
		}
		
	    }
	}
    }
    
#####################################

    my $refScoreH = undef; #hash-ref for storing selfScores
    $refScoreH = { ('dattype' => "selfScore") } if($selfcFl);
    
    my $HTMLstr = undef;
    my $alnData = undef; # for txt-output e.g. (structured array of alnData (only query lines)
    my $assignedHitsL = undef; # hash with key queryID, containing chainedBlastHits for assigned queryIDs
    my $exonFileD=undef; ## exonfile for Tracker program input ( listref of exon positions wrt chromosom))
    my $assignedRange=undef; ## sequence intervall where assigned hits lie
    my $ratios=undef; ### list of intra_inter_score_ratios for assigned hits
    ($HTMLstr, $alnData, $assignedHitsL, $exonFileD, $assignedRange, $ratios) = getHTMLstrRef($combPlocis, 
			     $cOptLociRef, 
			     $alnStrRef, 
			     $pfields,
			     $qfields,
			     # "syntenyRegion: <b>" . $blastFile ."</b>"
			     # . 
			     "<br/>AlnScore=<b>" . $alnScore . "</b>" . (($relScFl) ? " (based on relative scores)" : " (based on absolute scores)") 
			     . "<br/>AlnType: " . ( (${whichAlnInfo}) ? "reversed (decreasing chromosome positions)" : "normal (increasing chromosome positions)") . "\n"
			     . "<br/>QueryProtein(s)=<b>" . $refProtStr . "</b><br/>\n"
			     . "TargetGenome=" . $synOrg . "</b><br/>\n"
			     # . "<br/>qLoci-Score_CutOff: ${scutoff}\n"
			     # . "<br/>qLoci-MaxProtOverlap: ${protOvl}\n"
			     # . "<br/>loci-grouping type: " 
			     # . (($sFl) ? "mean cds/blastHit position" : "loci intervalls") . "<br/>\n"
			     , $regID  # for unique links in html after sorting
			     , $refOrgName # for links to Ensembl web site
			     , $ensemblQstr	### for ensembl link to query archive version
			     , $refScoreH # for getting selfcomparison-scores as hash for weighting of other blasthits
			     , 1 # $slicePos[0] # correction for getting real chromosome positions in tables instead of slice-positions (p-loci-proteins)
			     , $compAllFl # flag indicating to do the compara-test for every table entry
			     , $initF # optional initFile for EnsEMBL access methods
			     , $synOrg # needed targetName for comparaResult & Link
			     , $compData #contains all compara-homolog data with targOrg);
			     , $blHitsubDir
			     , $gRankData # contains global-best-hits with ranknumber for targOrg
			     , $refscores ## contains best and 2best scores of reference proteins against reference genome
			     , $refLidxs  
			     );
    
    my @PNGstrD = ();    # for a detailed png-graphic of alignment (with arrows and names)
    @PNGstrD = getPNGstrRef($combPlocis, 
			   $cOptLociRef, 
			   $alnStrRef,
			    $regID  # for unique links in html after sorting
			    , $whichAlnInfo,
			    700, # width in pixels
			    1, #detail flag
			    $dpColors, $alnIdxR
			    , $refLidxs);
    
    my @PNGstrS = ();  # the smaller one for overview   
    @PNGstrS = getPNGstrRef($combPlocis, 
			    $cOptLociRef, 
			    $alnStrRef,
			    $regID  # for unique links in html after sorting
			    , $whichAlnInfo,
			    100,
			    0,
			    $dpColors, $alnIdxR,  $refLidxs);
    
    my $pHTMLstr= undef; # additional table only for ploci (ppfields are displayed here)
    my $dummy=undef;
    ($pHTMLstr, $dummy) = getHTMLstrRef($combPlocis
    																, undef
    																, undef
    																, $ppfields
    																, undef
    																, "<h4>Data table for reference protein loci</h4>QueryProtein(s)=<b>" 
    																  . $refProtStr . "</b><br/>QueryGenome=<b>" . $refOrgName 
    																  . "</b><br/>\n" #(positions for chromosome/contig)<br/>
			       												, "00"
			       												, $refOrgName
			       												, $ensemblQstr
			       												, undef #, 1 #$slicePos[0]
															, undef
															, undef
															, undef
															, undef
															, undef
															, undef
															, undef
															, $refscores
															) if($regID == 1);
    
    push(@{$assignedRange}, $whichAlnInfo) if($assignedRange);
    
    
    
    my $RatioSum= 0;
    foreach (@{$ratios}) { $RatioSum += $_->[0]; }  ### calculate sum of ratios (for use as sorting criterium instead of alnScore) 

    my $focalRatioSum= 0;
    foreach (@{$ratios}) { $focalRatioSum +=  ( $_->[0] / (1. + $_->[1]) ); }  ### calculate weighted sum of ratios to focus on region around focal gene (for use as sorting criterium instead of alnScore) 

    
    @$ratios = ( [ (1e-2, 0 )] ) if(! @$ratios); ## set a negativ log_score value
    
    my $logRatioSum= 0;
    foreach (@{$ratios}) { $logRatioSum += log($_->[0]); }  ### calculate sum of log ratios -> product of ratios; (for use as sorting criterium instead of alnScore) 
     
   my $flogRatioSum= 0;
    foreach (@{$ratios}) { $flogRatioSum += ( log($_->[0]) / (1. + $_->[1]) ); }  ### calculate weighted sum of log ratios -> product of ratios; (for use as sorting criterium instead of alnScore) 
   #foreach (@{$ratios}) { $flogRatioSum += ( log($_->[0] / (1. + $_->[1]) ) ); }  ### calculate sum of log weighted ratios -> product of ratios; (for use as sorting criterium instead of alnScore) 
   
   # $avgRatio /= @{$ratios} if($ratios && @{$ratios});
    
    
   # my $avgRatio= 1;
    #foreach (@{$ratios}) { $avgRatio *= $_; }  ### calculate product of ratios (for use as sorting criterium instead of alnScore)
    
    my @retA = ($alnScore, undef, undef, undef, undef, undef, undef, undef, undef,undef, undef, undef, undef);
    $retA[1]= \@qloci if(@qloci);
    $retA[2]= $HTMLstr if($HTMLstr);
    $retA[3]= $alnData if($alnData);
    $retA[4]= [ (\@PNGstrS, \@PNGstrD) ] if(@PNGstrS);
    $retA[5]= $alnDPref if($alnDPref);
    $retA[6]= $pHTMLstr if($pHTMLstr);
    $retA[7]= $refScoreH if($refScoreH);
    $retA[8]= $assignedHitsL if($assignedHitsL);
    $retA[9]= $exonFileD if($exonFileD);
    $retA[10] = $assignedRange if($assignedRange);
    $retA[11]= $alnDPpref if($alnDPpref); ## printerFriendly DotPlot
    $retA[12]= $RatioSum if($RatioSum); 
    $retA[13]= int($focalRatioSum * 100 +.5)/100. if($focalRatioSum); 
    $retA[14]= int($logRatioSum * 100 +.5)/100.; # if($logRatioSum); 
    $retA[15]= int($flogRatioSum * 100 +.5)/100.; # if($flogRatioSum); 
    return(@retA);
}
#####################################################################################################################################################################


##########################################################################
# returns the png-data for the alignment graphics
##########################################################################
sub getPNGstrRef {
    my ($plo, $qlo, $alnstr, $regID, $revFl, $size, $detailFl, $dpColors, $alnIdxA,  $refLidxs) = @_;

    $size = 100 if(! defined($size));
    $detailFl=0 if(! defined($detailFl));
    my $maxRelSc = 1000;  ###

    my $WIDTH = $size;

### count number of gaps on the reference side (these elements can be drawed shorter)
    my $refGapNumber=0;
    for(my $li=0; $li< @{$alnstr}; ++$li) { # go through alnStr for a first look
	my $cP = $alnstr->[$li]->[0];
	my $cQ = $alnstr->[$li]->[1];
	my $cplo = undef;
	my $cplonumb = undef;
	if($cP ne "---") {
		$cplonumb = int(substr($cP, 1));
		$cplo = $plo->[$cplonumb -1];
  }
	my $cqlo =undef;
	$cqlo = $qlo->[int(substr($cQ, 1))-1] if($cQ ne "---");
#	++$refGapNumber if(($cqlo) && (! $cplo));
	++$refGapNumber if(! $cplo);
    }
    
    my $locSepSize = 2;
    my $locThickS=2;
    my $glocSize = 1 + $locSepSize;
    my $alnLength = @{$alnstr};
    my $locSize = int(($WIDTH -15 - ($refGapNumber * $glocSize) - ($locSepSize * ($alnLength +1))) / ($alnLength - $refGapNumber)) + 1;
    $locSize = 2*$locSepSize if($locSize < 2*$locSepSize);
    
    $WIDTH = ($alnLength + 1) * $locSepSize + ($alnLength - $refGapNumber) * $locSize + $refGapNumber * $glocSize +1 +16;
    my $yQloc = 4;
    my $yPloc = $yQloc + 8 + $locThickS;
    my @yDescr = (0,0,-9);
    @yDescr = ($yPloc + 8 + $locThickS + 4, $yPloc + 8 + $locThickS + 14, $yPloc + 8 + $locThickS -6 ) if($detailFl);

    
    my $HEIGHT = $yPloc + $yQloc + $yDescr[2] +9  + $locThickS + 1; #24;
    my $img = new GD::Image($WIDTH, $HEIGHT);
    
    # color allocation
    my $black   = $img->colorAllocate(0,0,0);  # first is background color
    my $white   = $img->colorAllocate(255,255,255);
    my $gray    = $img->colorAllocate(128,128,128);

    my $dred    = $img->colorAllocate(180,12,12);
    my $dred2     = $img->colorAllocate(255,100,0);  # 2nd color for hightlighting ref-protein

 
    my $red     = $img->colorAllocate(255,0,0);
    my $red2     = $img->colorAllocate(255,100,0);

    my $yellow  = $img->colorAllocate(255,255,0);

    my $dgreen  = $img->colorAllocate(49,87,78);
    my $green   = $img->colorAllocate(0,255,0);

    my $cyan    = $img->colorAllocate(0,255,255);
    my $blue    = $img->colorAllocate(0,0,255);
    
    
    my ($fivePos, $threePos)= (0,$WIDTH - 8);
    ($fivePos, $threePos)=($threePos,$fivePos) if($revFl);
    $img->string(gdTinyFont, $fivePos, $yQloc , "5'", $white);
    $img->string(gdTinyFont, $threePos, $yQloc , "3'", $white);

    $img->string(gdTinyFont, 0, $yPloc , "5'", $white);
    $img->string(gdTinyFont, $WIDTH - 8, $yPloc , "3'", $white);

    my $curXpos = $locSepSize+8;
    my $curYpos = $yQloc;

    my $mapAreaStr="<div><map name=\"regionID" . $regID . "map${detailFl}\">\n";
    
    my $secDl=0;
    my $countDownM= scalar(@{$alnIdxA});
    for(my $li=0; $li< @{$alnstr}; ++$li) { # go through alnStr again, now drawing the stuff
	   my $cP = $alnstr->[$li]->[0];
	   my $cQ = $alnstr->[$li]->[1];
	   my $cplo = undef;
	   my $isRefProt = undef;
	   if($cP ne "---") {
		  my $cplonumb = int(substr($cP, 1));
		  $cplo = $plo->[$cplonumb -1];
     
      if(defined($refLidxs) && (index("_" . join("_", @{$refLidxs}) . "_", "_${cplonumb}_") >= 0)) {
				### highlight reference locus (focal gene of interest)
	   	  $isRefProt=1;
      }
     }  
	   
	   
	   my $cqlo =undef;
	   $cqlo = $qlo->[int(substr($cQ, 1))-1] if($cQ ne "---");
	
	my $qcolor = $black;
	my $pcolor = $black;

	$qcolor = $dgreen if($cqlo);
	$pcolor = $dred if($cplo);
	$pcolor = $dred2 if($cplo && $isRefProt);
	
	my($ppi, $qqi, $mSc, $relSc)=(undef, undef, undef, undef);
	if(($cqlo) && ($cplo)) {
	    
	    ($ppi, $qqi, $mSc, $relSc) = @{$alnIdxA->[--$countDownM]};
	    die "ERROR: matching indices not in correct order!\n" unless((int(substr($cP, 1)) == $ppi)
									 && (int(substr($cQ, 1)) == $qqi) );
	if(defined($dpColors) && defined($dpColors->[0])) {
		    my $mm= ($maxRelSc / (@{$dpColors->[0]} -1));
		    if(defined($relSc)) {
			my $colIdx = int(($relSc -1)/$mm);
			$colIdx = @{$dpColors->[0]} -1 if($relSc > $maxRelSc);
			die "ERROR: ColorIndex out of Range! (getPNGstrRef, ${colIdx}, ${relSc}, ${maxRelSc})\n" unless(($colIdx< @{$dpColors->[0]}) && ($colIdx>=0));
			$green = $img->colorAllocate(@{$dpColors->[0]->[$colIdx]});
		    } else {
			$green = $white;
		    }
		}
	    
	    ($qcolor, $pcolor) = ($green, $red);
	    ($qcolor, $pcolor) = ($green, $red2) if($isRefProt);
	    
	}
	
	
	
	my $clocSize = $locSize;
	$clocSize = $glocSize if(($cqlo) && (! $cplo));




	$curYpos = $yQloc;
	my $method = 'filledRectangle';
	$img->$method($curXpos, $curYpos, $curXpos + $clocSize, $curYpos + $locThickS, $qcolor);
	
	my $pori53=0;
	my $pori35=0;
	my $qori53=0;
	my $qori35=0;
	if($cplo) {
	    my $cPdescr = $cP;
	    $cPdescr = $cQ . " - " . $cPdescr if($cqlo);
	    if($detailFl) {
		my @plkey = keys %{$cplo};
		(@plkey == 1) or die "Error - There is more than one Hash-Entry for a combined_P_locus (getPNGstrRef)!\n";
		my @pProtKeys = sort(keys %{$cplo->{$plkey[0]}});
		$cPdescr .= " (";
		my $cPshort = "";
		my $dummy=0;
		for(my $ii=0; $ii<@pProtKeys; $ii++) {
		    foreach my $cpe (@{$cplo->{$plkey[0]}->{$pProtKeys[$ii]}}) {
			$cPdescr .= "; " if($dummy);
			$cPshort .= ";" if(($dummy) && (length($cPshort)>1));
			++$dummy;
			$cPdescr .= $cpe->{'protID'};
			$cPdescr .= " [" . $cpe->{'geneName'} . "]" if(length($cpe->{'geneName'})> 1);
			$cPshort .= $cpe->{'geneName'} if(length($cpe->{'geneName'})> 1);

			++$pori53 if($cpe->{'ori'} > 0);
			++$pori35 if($cpe->{'ori'} < 0);
			if($cqlo) {
			    my $cQori= int(getQlocDat($cqlo, $cpe->{'protID'}, 'ori'));
			    ++$qori53 if($cQori > 0);
			    ++$qori35 if($cQori < 0);
			}
		    }
		}
		$cPdescr .= ")";
		if(defined($mSc)) {
		    $cPdescr .= " score=${mSc}";
		    $cPdescr .= ", rscore=${relSc}" if(defined($relSc));
		}

		$img->string(gdTinyFont, $curXpos + int(($clocSize - (length($cP)*5))/2), $yDescr[2] , $cP, $pcolor);	
		if(length($cPshort)>1) {
		    my $dcPshort= substr($cPshort, 0, SynBlast::MyUtils::min(length($cPshort), int(2*$clocSize / 5.)));
		    $img->string(gdTinyFont, $curXpos, $yDescr[$secDl] , $dcPshort, $pcolor);
		    if(($secDl) || ((length($cPshort)*5)> $clocSize)) { $secDl = ( $secDl ) ? 0 : 1; }
		}
		
	    }
	    
	    $mapAreaStr .= "<area shape=\"rect\" coords=\"" 
		. join(",", ($curXpos, $yQloc, $curXpos + $clocSize, $HEIGHT -1 )) #$curYpos + $locThickS))  
		. "\" href=\"#regionID_" . $regID . "_" . $cP . "\" alt=\"" . $cP . "\" title=\"" . $cPdescr . "\">\n";
	}
	
	$curYpos = $yPloc;
	$method = 'filledRectangle';
	$img->$method($curXpos, $curYpos, $curXpos + $clocSize, $curYpos + $locThickS, $pcolor);
	
	my $arrowSize=SynBlast::MyUtils::min(1+ int($locThickS/2), int($locSize/3));
	if($cplo && $detailFl) {  ### painting of the arrows for qloci
	    if($cqlo) {
		if((($qori53) && (! $revFl)) || (($qori35) && ($revFl))) {
		    $img->$method($curXpos + $clocSize - 2 * $arrowSize, $yQloc, $curXpos + $clocSize, $yQloc + $locThickS, $black);
		    my $poly = new GD::Polygon;
		    $poly->addPt($curXpos + $clocSize - 2 * $arrowSize, $yQloc +$arrowSize + $locThickS);
		    $poly->addPt($curXpos + $clocSize - 2 * $arrowSize, $yQloc -$arrowSize);
		    $poly->addPt($curXpos + $clocSize, $yQloc + int($locThickS/2));
		    my $arrcolor=$qcolor;
		    if(! $pori53) {
			$img->openPolygon($poly, $arrcolor);
		    } else {
			$img->filledPolygon($poly, $arrcolor);
		    }
		    
		}
		if((($qori35) && (! $revFl)) || (($qori53) && ($revFl))) {
		    $img->$method($curXpos, $yQloc, $curXpos + 2 * $arrowSize, $yQloc + $locThickS, $black);
		    my $poly = new GD::Polygon;
		    $poly->addPt($curXpos + 2 * $arrowSize, $yQloc +$arrowSize + $locThickS);
		    $poly->addPt($curXpos + 2 * $arrowSize, $yQloc -$arrowSize);
		    $poly->addPt($curXpos, $yQloc + int($locThickS/2));
		    my $arrcolor=$qcolor;
		    if(! $pori35) {
			$img->openPolygon($poly, $arrcolor);
		    } else {
			$img->filledPolygon($poly, $arrcolor);
		    }
		}
		if(! ($qori35 || $qori53)) { print STDERR "WARNING: There was a target-loci (qloci) without proper orientation!\n"; }
	    }


	    ###painting of the arrows for ploci
	    if($pori53) {
		$img->$method($curXpos + $clocSize - 2 * $arrowSize, $curYpos, $curXpos + $clocSize, $curYpos + $locThickS, $black);
		my $poly = new GD::Polygon;
		$poly->addPt($curXpos + $clocSize - 2 * $arrowSize, $curYpos +$arrowSize + $locThickS);
		$poly->addPt($curXpos + $clocSize - 2 * $arrowSize, $curYpos -$arrowSize);
		
		$poly->addPt($curXpos + $clocSize, $curYpos + int($locThickS/2));
		$img->filledPolygon($poly, $pcolor);
	    }

	    if($pori35) {
		$img->$method($curXpos, $curYpos, $curXpos + 2 * $arrowSize, $curYpos + $locThickS, $black);
		my $poly = new GD::Polygon;
		$poly->addPt($curXpos + 2 * $arrowSize, $curYpos +$arrowSize + $locThickS);
		$poly->addPt($curXpos + 2 * $arrowSize, $curYpos -$arrowSize);
		$poly->addPt($curXpos, $curYpos + int($locThickS/2));
		$img->filledPolygon($poly, $pcolor);
	    }
	    
	    if(! ($pori35 || $pori53)) { print STDERR "WARNING: There was a ref protein without proper orientation!\n"; }
	}
	

	
	$curXpos+=$clocSize + $locSepSize;
	
	
    }
    
    $mapAreaStr .= "</map></div>\n";
   
    my $retstr= $img->png;

    
    return(\$retstr, \$mapAreaStr);
}




##########################################################################
# returns the table row data for a certain p-locus assignment of alignment
##########################################################################
sub getPlocusHTML {
    my ($cplo, $cP, $cqlo, $cQ, $calnf, $alnDArrR, $plStyleStr, $qlStyleStr, $descrP, $descrQ, $regID, $refOrgName, $refScoreH, $sliceStart, $sliceAdapt, 
	$synOrg, $compData, $blHitHash, $blHitsubDir, $gRankData, $ensemblQstr, $exonL, $minmaxRange, $refscores, $ratios, $distToFocal) = @_; 
    
    
    my $HTMLs = "";
    my @cPalnData=();
    
    my $rowspanP = 1;
    @{$descrQ} = () if(! defined($cQ));
    
    my %ensemblRefTypes = ( 'protID' => 'protview?peptide=',
			    'geneID' => 'geneview?gene=',
			    'transID' => 'transview?transcript='
			    );
    my $ensemblSite = "http://www.ensembl.org";   # should correspond to archive-version of info-file version
    $ensemblSite = "http://" . $ensemblQstr . ".archive.ensembl.org" if(defined($ensemblQstr));
    # "http://aug2006.archive.ensembl.org"
    
    
    my @Pcells=();
    my @Qcells=();
    
    
    if(defined($cplo) && @{$descrP}) {
      
	my @plkey = keys %{$cplo};
	(@plkey == 1) or die "Error - There is more than one Hash-Entry for a combined_P_locus (getPlocusHTML)!\n";
	my @pProtKeys = sort(keys %{$cplo->{$plkey[0]}});
      $rowspanP = @pProtKeys;
      
    
      for(my $ii=0; $ii<@pProtKeys; $ii++) {
	my @cPalnSubDat=();
	my @tmpR=();
	my @tmpQ=();
	my $scounter=0;
	push(@cPalnSubDat, $cP);
	push(@cPalnSubDat, $cQ);

	
	foreach my $cpe (@{$cplo->{$plkey[0]}->{$pProtKeys[$ii]}}) {
	  ++$scounter;
	  if($scounter == 1) {
	    if(defined($cpe->{'protID'})) {
	      push(@cPalnSubDat, $cpe->{'protID'});
	    } else {
	      push(@cPalnSubDat, "unknownQueryID");
	    }
	  }
	  for(my $di=0; $di<@{$descrP}; $di++) {
	   
		if($descrP->[$di] eq "selfscore") {
			my $sesc=undef;
			$sesc= $refscores->{$cpe->{'protID'}}->[0] if(defined($cpe->{'protID'}) && defined($refscores->{$cpe->{'protID'}}));
			my $seva="not used";
			$seva=$sesc if(defined($sesc));
			$tmpR[$di] .= $seva;
			push(@cPalnSubDat, $seva) if($scounter == 1);
			$tmpR[$di] .= "<br/>\n" if(@{$cplo->{$plkey[0]}->{$pProtKeys[$ii]}} > 1);
		} elsif(defined($cpe->{$descrP->[$di]})) {
				push(@cPalnSubDat, $cpe->{$descrP->[$di]}) if($scounter == 1);
				if(exists($ensemblRefTypes{$descrP->[$di]})) {
					my $curDat=$cpe->{$descrP->[$di]};
					$tmpR[$di] .= "<a href=\"" 
						. "${ensemblSite}/${refOrgName}/" 
						. $ensemblRefTypes{$descrP->[$di]}
						. $curDat
						. "\" target=\"_ensembl\">"; 
					$tmpR[$di] .= $cpe->{$descrP->[$di]};
					$tmpR[$di] .= "</a>\n";
				} else {
					my $value=$cpe->{$descrP->[$di]};
					if(index($descrP->[$di], "Pos", 1)>=0) {
						# correct positions to correspond to real chromosome positions
						$value += $sliceStart -1;
					}
					$tmpR[$di] .= $value; #$cpe->{$descrP->[$di]};
				}
				$tmpR[$di] .= "<br/>\n" if(@{$cplo->{$plkey[0]}->{$pProtKeys[$ii]}} > 1);
		} else {
			$tmpR[$di] .= "noData";
			push(@cPalnSubDat, "noData") if($scounter == 1);
		}
	  }

	  if(defined($cqlo)) { ## exonfile
	    my $pName=$cpe->{'protID'};
	    $pName=$cpe->{'protName'} if($cpe->{'protName'} ne "-");
	    my $gName=$cpe->{'geneID'};
	    $gName=$cpe->{'geneName'} if($cpe->{'geneName'} ne "-");
	    my $tHits=getQlocDat($cqlo, $cpe->{'protID'}, "hits");
	    if($tHits ne "NoData") {
	      my @curExon=();
	      push(@curExon, "${pName} ${gName}");
	      push(@curExon, $cpe->{'protID'});
	      @{$minmaxRange} = (SynBlast::MyUtils::min($minmaxRange->[0], getQlocDat($cqlo, $cpe->{'protID'}, "startPos"))
				 , SynBlast::MyUtils::max($minmaxRange->[1], getQlocDat($cqlo, $cpe->{'protID'}, "endPos")));
	      
	      foreach my $thi (@{$tHits}) {
		push(@curExon, [ ($thi->[8], $thi->[9]) ]);
	      }
	      push(@{$exonL}, \@curExon);
	    }
	  }
	  	
	  if(defined($cqlo) && @{$descrQ}) {
	    for(my $di=0; $di<@{$descrQ}; $di++) {
	      
	      
	      ################ h i t s ############################
	      if($descrQ->[$di] eq "hits") {
			  my $tHits=getQlocDat($cqlo, $cpe->{'protID'}, $descrQ->[$di]);
			  my $hitCount=0;
			  if($tHits eq "NoData") {
			    $tmpQ[$di] .= $tHits;
			  } else {
			    $hitCount=@{$tHits};
			    if(defined($blHitHash)){
			      $tmpQ[$di] .= "<a href=\"./${blHitsubDir}/rID" . $regID . "_" . $cpe->{'protID'} . "." . $synOrg . ".chainedBlastHits.txt" . "\" target=\"_hits_" . $cpe->{'protID'} . "\">";
			      $tmpQ[$di] .= "${hitCount}";
			      $tmpQ[$di] .= "</a>\n";				
			      if(exists($blHitHash->{$cpe->{'protID'}})) {
				warn "WARNING: blastHits for this queryID will be overwritten!\n";
			      }
			      $blHitHash->{$cpe->{'protID'}} = $tHits;		
			      
			    } else {
			      $tmpQ[$di] .= "<table><tbody>\n";
			      foreach my $thi (@{$tHits}) {
				$tmpQ[$di] .= "<tr><td>" . join("</td><td>", @{$thi}) . "</td></tr>\n"
				  }
			      $tmpQ[$di] .= "</tbody></table>\n";
			    }
			  }
			  push(@cPalnSubDat, "${hitCount}") if($scounter == 1);
	   }   ################ h i t s  end ############################
	   
	   
	   
	   ################ c o m p a r a ############################
	      elsif($descrQ->[$di] eq "compara") {
		
			    my $coordType = ""; #"scaffold"; #empty string works for all (other) contig/scaffold/group types...
			    my $ccontig=getQlocDat($cqlo, $cpe->{'protID'}, 'chr');
			    if(defined($sliceAdapt) && ($ccontig ne "NoData")) { # && ($scounter == 1)) {
						my $cstartp=getQlocDat($cqlo, $cpe->{'protID'}, 'startPos');
						my $cendp=getQlocDat($cqlo, $cpe->{'protID'}, 'endPos');
						my $cRefGID=$cpe->{'geneID'};
						
						### it seems not to work properly in all cases, so leave coordType undefined/empty and EnsemblAPI automatically take toplevel type
						#if(int($ccontig)) {
							#$coordType = "chromosome";
						#}
						my ($txtResult, $resHash, $gid2Name) = SynBlast::EnsEMBLaccess::checkComparaData($compData->{$cRefGID}, $sliceAdapt, $synOrg, $coordType, $ccontig, $cstartp, $cendp, $cRefGID, 0, 1);
				
						chomp($txtResult);
						my $compResHtml = SynBlast::MyUtils::comparaRes2HTML($resHash, $synOrg, $ensemblSite, $gid2Name);
						$tmpQ[$di] .= $compResHtml;
						$tmpQ[$di] .= "<br/>\n" if(@{$cplo->{$plkey[0]}->{$pProtKeys[$ii]}} > 1);
						push(@cPalnSubDat, $txtResult) if($scounter == 1);
			      } else {
							$tmpQ[$di] .= "-";
							$tmpQ[$di] .= "<br/>\n" if(@{$cplo->{$plkey[0]}->{$pProtKeys[$ii]}} > 1);
							push(@cPalnSubDat, "noData") if($scounter == 1);
			      }
			    
			  } ################ c o m p a r a  end ############################
			 
			   elsif($descrQ->[$di] eq "intra_inter") {
				my $intrainterR="noData";
				 if(exists($refscores->{$cpe->{'protID'}})) {
					 my @refscores= @{$refscores->{$cpe->{'protID'}}};
					 my ($intraref, $inter, $quot);
					 $refscores[1] = 0 if(! $refscores[1]);
					 #$refscores[1] = $refscores[0] - 0.01 if($refscores[0] && ($refscores[0] - $refscores[1] <= 0));  ### limit intrascore by 1percent
					 $intraref =  int(($refscores[0] - $refscores[1]) / $refscores[0] * 100)/100. if($refscores[0] && defined($refscores[1]));
					 #$intraref = 0.01 if((! $intraref) || ($intraref <= 0));
					 $intraref = 0.01 if(defined($intraref) && ($intraref <= 0));
					 
					 my $csco= getQlocDat($cqlo, $cpe->{'protID'}, 'score');
					  $inter = int(($refscores[0] - $csco) / $refscores[0] * 100)/100. if($refscores[0] && ($csco ne "NoData"));
					  if(defined($inter)) { ## only if data was available
						$inter = 0.01 if((! $inter) || ($inter <= 0));    # if target score was larger equal bestRefScore -> set inter to zero, too  ## but limit by 1 percent
						$intrainterR= "${intraref} / ${inter} = ";
					  	#if(!($inter) && $intraref) {
						#	$quot= 1000;
						  #} elsif
						  if($inter && defined($intraref) ) {
							  $quot = int(100. * $intraref / $inter +.5)/100.;
							  $quot = 0.01 if((! $quot) || ($quot < 0.01));  ## limit ratio by 1 percent 
						  }
					  }
					  if(defined($quot)) {
						  $intrainterR .= $quot;
						  push(@{$ratios}, [ ($quot, $distToFocal) ] );
					  }
					  
				} 
				 $tmpQ[$di] .= $intrainterR;
			        $tmpQ[$di] .= "<br/>\n" if(@{$cplo->{$plkey[0]}->{$pProtKeys[$ii]}} > 1);
			     push(@cPalnSubDat, $intrainterR) if($scounter == 1);
			    
			  } ################ intra_inter end ############################
			 
			 
			   elsif($descrQ->[$di] eq "globalRank") {
				   my $glRankRes="";
				   my $tmpQchr=getQlocDat($cqlo, $cpe->{'protID'}, 'chr');				   
				   if(exists($gRankData->{$cpe->{'protID'}}->{$tmpQchr})) {
				   		#### for current assignments region is a list of nbest global hits available...
				   		# go through list and check on overlap in positions (and identical score)	
				   		my $cListRef=$gRankData->{$cpe->{'protID'}}->{$tmpQchr};
				   		my $hitfoundFl=0;
				   		my $globHitIdx=0;
				   		my $tmpQstart=getQlocDat($cqlo, $cpe->{'protID'}, 'startPos');				   
				   		my $tmpQend=getQlocDat($cqlo, $cpe->{'protID'}, 'endPos');				   
				   		while((! $hitfoundFl) && ($globHitIdx < @{$cListRef})) {
				   			my $cHiStart=$cListRef->[$globHitIdx]->[1];
				   			my $cHiEnd=$cListRef->[$globHitIdx]->[2];
				   			my $cHiScore=$cListRef->[$globHitIdx]->[4];
				   			if(SynBlast::MyUtils::getOverlapSize($cHiStart, $cHiEnd, $tmpQstart, $tmpQend)) {
				   				#### hit found due to overlap?;
				   				#$glRankRes .= "" . ($cListRef->[$globHitIdx]->[0]);
				   				++$hitfoundFl;
									if(int($cHiScore) == int(getQlocDat($cqlo, $cpe->{'protID'}, 'score'))) {
										### score identical too?
										$glRankRes .= "" . ($cListRef->[$globHitIdx]->[0]); # rank
									} else {
											$glRankRes .= "(" . ($cListRef->[$globHitIdx]->[0]) . ")"; # rank in brackets
									}
								}
				   			++$globHitIdx;
				   		} # while
				   		$glRankRes .= "noData" if(! $hitfoundFl);
				   	} else {
				   		$glRankRes .= "noData";
				   	}
				   $tmpQ[$di] .= $glRankRes;
			     $tmpQ[$di] .= "<br/>\n" if(@{$cplo->{$plkey[0]}->{$pProtKeys[$ii]}} > 1);
			     push(@cPalnSubDat, $glRankRes) if($scounter == 1);
			}			
			 ############# normal q-locus entry #####################
			else {
			  my $tmpQD=getQlocDat($cqlo, $cpe->{'protID'}, $descrQ->[$di]);
			  $tmpQ[$di] .= $tmpQD;
			  $tmpQ[$di] .= "<br/>\n" if(@{$cplo->{$plkey[0]}->{$pProtKeys[$ii]}} > 1);
			  push(@cPalnSubDat, $tmpQD) if($scounter == 1);
			}
	    }
		    
		    
		    
		    
	    ###
	    $refScoreH->{$cpe->{'protID'}} = getQlocDat($cqlo, $cpe->{'protID'}, 'score') if(defined($refScoreH));
	  } else {
	    $refScoreH->{$cpe->{'protID'}} = 0 if(defined($refScoreH));
	    ###
	  }
		
	}
	push(@Pcells, \@tmpR);
	push(@Qcells, \@tmpQ) if(defined($cqlo) && @{$descrQ});
	push(@cPalnData, \@cPalnSubDat);
      }
	
    }
    
    $HTMLs .= "<tr>";
    $HTMLs .= "<td rowspan=\"" . $rowspanP . "\"${plStyleStr}><a name=\"regionID_" . $regID . "_" . $cP . "\"";
    $HTMLs .= " href=\"#regionID_00_${cP}\"" if(($cP ne "---") && ($regID ne "00"));
    $HTMLs .= ">" . $cP . "</a></td>\n";
    
    if(defined($cQ)) {
	$HTMLs .= "<td rowspan=\"" . $rowspanP . "\"${qlStyleStr}>";
	$HTMLs .= "<a name=\"regionID_" . $regID . "_" . $cQ . "\"></a>" if(! scalar(@{$descrQ}));
	$HTMLs .= $cQ . "</td>\n";
	$HTMLs .= "<td rowspan=\"" . $rowspanP . "\">" . $calnf . "</td>\n";
    }
    
    for(my $lu=0; $lu<@{$descrP}; $lu++) {
	my $colstyle="";
	my $cDat=((@Pcells) ? $Pcells[0]->[$lu] : "-");
	if(($descrP->[$lu] eq "optQloci") && ($cDat ne $cQ) && ($cDat ne "-") && ($cDat ne "noMatch")) {
	  $colstyle=" class=\"notOptLoci\""; # style=\"background-color: rgb(128, 80, 80);vertical-align: center;\"";
	}
	$HTMLs .= "<td${colstyle}>" . $cDat . "</td>\n";
	
    }
    
    
    $HTMLs .= "<td rowspan=\"" . $rowspanP . "\"${qlStyleStr}><a name=\"regionID_" . $regID . "_" . $cQ . "\"></a>" . $cQ . "</td>\n" if(defined($cQ) && scalar(@{$descrQ}));
    for(my $lu=0; $lu<@{$descrQ}; $lu++) {
	$HTMLs .= "<td>" . ((@Qcells) ? $Qcells[0]->[$lu] : "-") . "</td>\n";
    }
    
    $HTMLs .= "</tr>";
    
    for(my $la=1; $la<$rowspanP; $la++) {
	$HTMLs .= "<tr>";
	for(my $lu=0; $lu<@{$descrP}; $lu++) {
	    my $colstyle="";
	    my $cDat= $Pcells[$la]->[$lu];
	    if(($descrP->[$lu] eq "optQloci") && ($cDat ne $cQ) && ($cDat ne "-") && ($cDat ne "noMatch")) {
	      $colstyle=" class=\"notOptLoci\""; #" style=\"background-color: rgb(128, 80, 80);vertical-align: center;\"";
	    } 
	    $HTMLs .= "<td${colstyle}>" . $cDat . "</td>\n";
	  }
	
	for(my $lu=0; $lu<@{$descrQ}; $lu++) {
	    $HTMLs .= "<td>" . ((@Qcells > $la) ? ($Qcells[$la]->[$lu]) : "-") . "</td>\n";
	}
	$HTMLs .= "</tr>";
      }
    


    push(@{$alnDArrR}, \@cPalnData) if(defined($alnDArrR));
    return(\$HTMLs);
  }    






##########################################################################
# returns an alignment result block including HTML-assignment-table as string ref containing/describing the 
# selected alignment including the specified table-columns for p- and q-loci
##########################################################################
sub getHTMLstrRef {
    my ($plo, $qlo, $alnstr, $descrP, $descrQ, $tabDescr
		, $regID, $refOrgName, $ensemblQstr, $refScoreH, $sliceStart
		, $compAllFl, $initF, $targetOrg, $compData, $blHitSubDir, $gRankData, $refscores, $refLidxs) = @_; 

    $compAllFl=0 if(! defined($compAllFl));
    
    
    # is used to print the real position on the chromosome of the cds region (if 1 or not defined, the slice positions are displayed)
    $sliceStart=1 if(! defined($sliceStart)); 

    
    @{$descrP} = qw(protID geneName optQloci) if(! defined($descrP));
    @{$descrQ} = qw() if(! defined($descrQ));
    
    my @descrPn;
     @descrPn= @{$descrP};
    push(@descrPn, "selfscore") if(defined($refscores) && (! defined($alnstr)));  ## only print selfscore in refloci table
      
    
    my $plStyleStr=""; #" style=\"background-color: rgb(255, 255, 102);vertical-align: center;\""; ## is set now via css-class property and css.style in html.css
    my $qlStyleStr=""; #" style=\"background-color: rgb(51, 255, 51);vertical-align: center;\"";
    
    my $HTMLs = ""; #"<div class=\"aliresult\">\n";

    my @alnDArr=();
    my %blastHits=();
    
    my @ratios=();
    
    my @exonL=(); 
    my @minmaxRange=(undef, undef);

    $HTMLs .= ($tabDescr . "\n") if(defined($tabDescr));
    $HTMLs .= "<table class=\"aliresult\">\n<tr>";
    $HTMLs .= "<th class=\"refloci\"${plStyleStr}>ref-loci</th>\n" if(defined($plo));
    $HTMLs .= "<th class=\"targetloci\"${qlStyleStr}>target-loci</th>\n" if(defined($qlo));
    $HTMLs .= "<th class=\"targetloci\"${qlStyleStr}>aln</th>\n" if(defined($alnstr));
    $HTMLs .= "<th class=\"refloci\"${plStyleStr}>" . join("</th><th class=\"refloci\"${plStyleStr}>", @descrPn) . "</th>\n" if(@descrPn);
    $HTMLs .= "<th class=\"targetloci\"${qlStyleStr}>target-loci</th>\n" if(@{$descrQ});
    $HTMLs .= "<th class=\"targetloci\"${qlStyleStr}>" . join("</th><th class=\"targetloci\"${qlStyleStr}>", @{$descrQ}) . "</th>\n" if(@{$descrQ});
    $HTMLs .= "</tr>\n";
    
    
    if(defined($alnstr)) {
	my $slice_adaptor = undef;
	if($compAllFl) {
	    #get connection to core-DB of current organism...
	    $slice_adaptor = SynBlast::EnsEMBLaccess::getDBadaptor($targetOrg, "core", "Slice");
	}
	
	
	
	for(my $li=0; $li< scalar(@{$alnstr}); ++$li) { # go through alnStr
	    my $cP = $alnstr->[$li]->[0]; # e.g. p1 or ---
	    my $cQ = $alnstr->[$li]->[1]; # e.g. q4 or ---
	    my $cplo = undef;
	    $cplo = $plo->[int(substr($cP, 1))-1] if($cP ne "---");    #ploci-data
	    my $cqlo =undef;
	    $cqlo = $qlo->[int(substr($cQ, 1))-1] if($cQ ne "---");  # qloci-data
	   my $distToFocal=undef;
	   if($cplo) {
		   $distToFocal= abs(int(substr($cP, 1)) - int($refLidxs->[0]));
	   }
	   #print "focalGeneList=" . join(", ", @$refLidxs) . "\tDistance to Focal=" . $distToFocal . "\n"; # unless($distToFocal);
	   
	    $HTMLs .= ${getPlocusHTML($cplo, $cP, $cqlo, $cQ, $alnstr->[$li]->[2], \@alnDArr, " class=\"refloci\"${plStyleStr}", " class=\"targetloci\"${qlStyleStr}", \@descrPn, $descrQ, $regID, $refOrgName, $refScoreH, $sliceStart, $slice_adaptor, $targetOrg,$compData, \%blastHits, $blHitSubDir, $gRankData, $ensemblQstr, \@exonL, \@minmaxRange, $refscores, \@ratios, $distToFocal)};
	    
	}
      } elsif(defined($plo)) {
	for(my $li=0; $li< scalar(@{$plo}); ++$li) { # go through ploci
	  $HTMLs .= ${getPlocusHTML($plo->[$li]
					, "p" . ($li + 1)
					, undef, undef, undef, undef
					, " class=\"refloci\"${plStyleStr}"
					, undef
					, \@descrPn
					, undef
					, $regID
					, $refOrgName
					, undef
					, $sliceStart
					, undef
					, undef
					, undef
					, $ensemblQstr
					, undef
					, undef
					, $ensemblQstr
					, undef
					, undef
					, $refscores
					, undef
					, undef
					)};
	}
      } elsif(defined($qlo)) {
	for(my $li=0; $li< scalar(@{$qlo}); ++$li) { # go through qloci
	  $HTMLs .= ${getPlocusHTML(undef, undef, $qlo->[$li], "q" . ($li + 1), undef, undef, undef, " class=\"targetloci\"${qlStyleStr}", undef, $descrQ, $regID, $refOrgName, undef, $sliceStart, undef, undef, undef, undef, undef, undef)};
	}
    } else {
      die("\nNo valid parameters used! (getHTMLstrRef)\n");
    } 

    $HTMLs .= "</table>\n"; #</div>\n";

    return(\$HTMLs, \@alnDArr, \%blastHits, \@exonL, \@minmaxRange, \@ratios);
}




##########################################################################
# returns the q-locus assignment data for a certain reference locus
##########################################################################
sub getQlocDat {
    my ($qloc, $queryID, $ekey) = @_;
    my @qlkey = keys %{$qloc};
    (@qlkey == 1) or die "Error - There is more than one Hash-Entry for a combined_Q_locus (getQlocDat)!\n";
    if(defined($qloc->{$qlkey[0]}->{$queryID}->[0]->{$ekey})) { # hitdata for queryID existing...
	return($qloc->{$qlkey[0]}->{$queryID}->[0]->{$ekey});
    }
    return("NoData");
}


##########################################################################
# returns the p-locus data 
##########################################################################
sub getPlocDat {
    my ($ploc, $entryNumber, $ekey) = @_;
    my @plkey = keys %{$ploc};
    (@plkey == 1) or die "Error - There is more than one Hash-Entry for a combined_P_locus (getPlocDat)!\n";
    my $returnStr="";
    my @skeys = (sort keys %{$ploc->{$plkey[0]}});
    if(($entryNumber < @skeys) && defined($ploc->{$plkey[0]}->{$skeys[$entryNumber]}->[0]->{$ekey})) { # ploci - data for queryID existing...
	return($ploc->{$plkey[0]}->{$skeys[$entryNumber]}->[0]->{$ekey});
    }
    #DataBrowser::browse($ploc->{$plkey[0]});
    return("NoData");
}





##########################################################################
#  get list of optimal chainedTargetLoci
##########################################################################
sub getOptHitList {
    my ($refexons, $locisref) = @_; 
    my @outp=();
    for(my $li=0; $li<@{$refexons}; $li++) {

#	push(@outp, [ ( "p" . ($li+1) . "=" . $refexons->[$li]->[2] , getOptLociStr($refexons->[$li]->[2], $locisref) ) ] );

	my @cListList=();
	my @plkey = keys %{$refexons->[$li]};
	(@plkey == 1) or die "Error - There is more than one Hash-Entry for a combined_P_locus (getOptHitList)!\n";
	my $cplociHash=$refexons->[$li]->{$plkey[0]};
	foreach my $cplociE (keys %{$cplociHash}) {
	    foreach my $cPE (@{$cplociHash->{$cplociE}}) {
		my $protID=$cPE->{'protID'};
		my $protName=$cPE->{'protName'};
		$protName= ($protName ne "-") ? ($protID . "=>" . $protName) : $protID;
		my ($optScoreV, $qsRef) = getOptLociDat($protID, $locisref);
		$cPE->{'optQloci'} = ((@{$qsRef}) ? join(",", @{$qsRef}) : "noMatch");
		
		push(@cListList, [ ( "p" . ($li+1) . "=>" . $protName , ("maxScore=" . $optScoreV . "(" . join(",", @{$qsRef}) . ")") ) ] );
		#	    push(@outp, [ ( "p" . ($li+1) . "=>" . $protID , getOptLociStr($protID, $locisref) ) ] );
	    }
	}
	push(@outp, \@cListList);
	
    }
    (@outp == @{$refexons}) or die "Error - the optHitList should be of the same length as ploci-number! (getOptHitList)\n";
    return(\@outp);
}

sub getOptLociDat {
    my ($exonID, $locisref) = @_; 
    my $curOptV=0;
    my @curOptLociN=();
    for(my $ll=0; $ll<@{$locisref}; $ll++) {
	my @lkey= keys %{$locisref->[$ll]};
	(@lkey == 1) or die "Error: There are more than one Hash-Entries for a putative loci (getOptLociStr)\n";
	if(exists($locisref->[$ll]->{$lkey[0]}->{$exonID})) {
	    # match in current loci?
	    #my $curV=$locisref->[$ll]->{$lkey[0]}->{$exonID}->[-1];
	    die "ERROR: combined_q_Loci should have not more than one Entry for the same (combined) range and ProteinQueryID! (" 
		. @{$locisref->[$ll]->{$lkey[0]}->{$exonID}} . ")\n" unless ( @{$locisref->[$ll]->{$lkey[0]}->{$exonID}} == 1);
	    foreach my $cE (@{$locisref->[$ll]->{$lkey[0]}->{$exonID}}) {
		my $curV=$cE->{'score'};
		if($curV == $curOptV) {
		    push(@curOptLociN, $ll);
		} else {
		    if($curV > $curOptV) {
			$curOptV=$curV;
			@curOptLociN=( $ll );
		    }
		}
	    }
	}
    }
    if(@curOptLociN) {
	my @outf= ($curOptV);
	my @qs=();
	for(my $ll=0; $ll<@curOptLociN; $ll++) {
	    push(@qs, "q" . ($curOptLociN[$ll] +1));
	    #	    $outstr= $outstr . ", " if($ll<$#curOptLociN);
	}
	push(@outf, \@qs);
	return(@outf);
    }
    return(("noMatch", [ () ]));
}









##########################################################################
# alignLocis - main sub for the Needleman-Wunsch-like gene order alignment, 
# includes generation of dotplot-graphics png-data (if flag doDPFl set)
##########################################################################
sub alignLocis {
  my ($refexons
      , $locisref
      , $gapsimP
      , $gapsimQ
      , $mismV
      , $verboseFl
      , $penaltyFl
      , $penaltyFa
      , $regID
      , $dpColors
      , $printcFl
      , $myScoreFunction
      , $relAlScFl
      , $refLidxs
      , $doDPFl
      , $maxDist
      , $blqueryStI
      , $bltargetStI
      , $regDescr
      , $refDescr
      , $refscores
      ) = @_;
  
    $verboseFl=0 if(! defined($verboseFl));
    $penaltyFl=0 if(! defined($penaltyFl));
    $penaltyFa=10 if(! defined($penaltyFa));
    $relAlScFl=0 if(!defined($relAlScFl)); # flag for use of relativeScore as matchSimilarity instead of absolute score
    $doDPFl = 0 if(! defined($doDPFl)); # do not create dotplot data per default
    $myScoreFunction = \&getAlnLociScore if(! defined($myScoreFunction));

    $regDescr="target" if(! defined($regDescr));
    $refDescr="ref" if(! defined($regDescr));
    
    my $maxRelSc=1000; #used to scale the relative score and color allocation
       
   my @colDistances=();  ## to keep the distances between the target loci columns
    for(my $ii=1; $ii< @{$locisref}; ++$ii) {
      my @lkey1= keys %{$locisref->[$ii]};
      (@lkey1 == 1) or die "Error: There are more than one Hash-Entries for a putative locus\n";
      my @lkey0= keys %{$locisref->[$ii-1]};
      (@lkey0 == 1) or die "Error: There are more than one Hash-Entries for a putative locus\n";
      my($k1from,$k1to)=split(/_/, $lkey1[0]);
      my($k0from,$k0to)=split(/_/, $lkey0[0]);
      my $kdist= $k1from - $k0to;
      $kdist= $k0from - $k1to if($k0from > $k1to);
      push(@colDistances, $kdist);
      #print "Distance col ${ii} to col " . (${ii}-1) . " = ${kdist} (keys were " . $lkey1[0] . " " . $lkey0[0] . ")\n";
    }
  
  my $maxColDis=1;
  $maxColDis=SynBlast::MyUtils::max(@colDistances) - $maxDist if(defined($maxDist) && @colDistances);
  my @colDistSepSizes=();
  my $maxDistSepS=undef;
  
  
  ###  ### simple dotplot png-str
    # image setup
    my($locSepSize, $locSize, $imgDP, $black, $lblack, $gray, $lgray, $dred, $red, $dgreen, $green, $yellow, $white, @myfont, @myfontb,  @myfonth, $eehei); 
    my($imgDPp, $pblack, $plblack, $pgray, $plgray, $pdred, $pred, $pdgreen, $pgreen, $pyellow, $pwhite, $ehei); 
    my $descrFl=10;  ### additionally write query gene names in rows (of max this length (characters))
    my $descrL=0;

  if($doDPFl) {
      $locSepSize = 2;
      my $qlength= @{$locisref};
      my $plength= @{$refexons};
      $locSize = 8;
      $locSize = 16 if($printcFl);
      $locSize = 2*$locSepSize if($locSize < 2*$locSepSize);
      my $myfonth=$locSize+1;
      my ($myfont, $myfontb);
      my @myfonts=(gdGiantFont, gdLargeFont, gdMediumBoldFont, gdSmallFont, gdTinyFont);
      @myfonts= reverse @myfonts;
      while(@myfonts) {
	      $myfont=shift(@myfonts);
	      $myfontb=$myfont->width;
	      $myfonth=$myfont->height;
		if( ($myfont == gdTinyFont) || (($myfonth<= $locSize) && ($myfontb <= $locSize)) ) {
			push(@myfont, $myfont);
			push(@myfontb, $myfontb);
			push(@myfonth, $myfonth);
		}
      }
      
      
      
      my $WIDTH = $qlength * ($locSize + $locSepSize) + $locSepSize +1;
      
      
      $descrL=$descrFl * $myfontb[-1] + 2* $locSepSize if($descrFl);
      
      $WIDTH +=  $descrL;
      
      my $HEIGHT = $plength * ($locSize + $locSepSize) + $locSepSize +1;

	$eehei=3* $locSize + 2*$locSepSize;
	$imgDP = new GD::Image($WIDTH, $HEIGHT + $eehei);
       #$imgDP = new GD::Image($WIDTH, $HEIGHT);
      # color allocation
      $black   = $imgDP->colorAllocate(0,0,0);  ## first is background color
      $lblack   = $imgDP->colorAllocate(62,62,62);
      $white   = $imgDP->colorAllocate(255,255,255);
      $gray    = $imgDP->colorAllocate(128,128,128);
      $lgray    = $imgDP->colorAllocate(168,168,168);
      $dred    = $imgDP->colorAllocate(180,12,12);
      $red     = $imgDP->colorAllocate(255,0,0);
      $yellow  = $imgDP->colorAllocate(255,128,0);
    	$dgreen  = $imgDP->colorAllocate(12,180,12);
      $green   = $imgDP->colorAllocate(0,255,0);
   #my $cyan    = $imgDP->colorAllocate(0,255,255);
    #my $blue    = $imgDP->colorAllocate(0,0,255);
######
  if($printcFl) {
	
	my $addWidth=0;
	
	if(defined($maxDist)) {
		$maxDistSepS=int($locSize/2.)+1;
		my $ffac= ($maxDistSepS)/ $maxColDis / $maxColDis;
		for(my $dli=0; $dli< @colDistances; ++$dli) {
			my $curSepThick= 0;
			$curSepThick= 1+ int($ffac * ($colDistances[$dli] - $maxDist) * ($colDistances[$dli] - $maxDist)) if($colDistances[$dli] > $maxDist);
			push(@colDistSepSizes, $curSepThick);
			$addWidth += $curSepThick;
		}
	}
  
	#$ehei=$locSize + 2*$locSepSize;
	 $ehei=3* $locSize + 2*$locSepSize;
  $imgDPp = new GD::Image($WIDTH + $addWidth, $HEIGHT + $ehei);
	
  
	### color allocation, first is background color
	$plgray    = $imgDPp->colorAllocate(168,168,168);
	$pblack   = $imgDPp->colorAllocate(0,0,0);
	$pwhite   = $imgDPp->colorAllocate(255,255,255);
	#$plblack   = $imgDPp->colorAllocate(62,62,62);
	$pgray    = $imgDPp->colorAllocate(128,128,128);
	$pdred    = $imgDPp->colorAllocate(180,12,12);  ### color for printerFriendlyFl
	   $pred     = $imgDPp->colorAllocate(255,0,0);
      $pyellow  = $imgDPp->colorAllocate(255,128,0);
	$pdgreen  = $imgDPp->colorAllocate(12,180,12);
	$pgreen   = $imgDPp->colorAllocate(0,255,0);

# make the background transparent and interlaced
    #$imgDPp->transparent($white);
    #$imgDPp->interlaced('true');
      	# Put a black frame around the picture
        #$imgDPp->rectangle(0,0,$WIDTH-1,$HEIGHT-1,$pblack);
	
	
	
	
  }

    }
    
    my @matrix=();
    #initialize alignment_matrix...
    ## with zero-values (for endfree-gaps)...
    for(my $ii=0; $ii <= @{$refexons}; $ii++) {
	my @curLi=();
	push(@curLi, 0); #$ii * $gapsim); #0);
	if ($ii == 0) { 
	    for(my $jj=0; $jj < @{$locisref}; $jj++) {
		push(@curLi, 0); #($jj +1) * $gapsim);
	    }
	} 
	push(@matrix, \@curLi);
    }

    # init score_matrix ##for better performance while backtracking etc., contains matchSims
    my @scmatrix=();
    
    for(my $ii=1; $ii <= @{$refexons}; $ii++) {
	my @curLi=();
	for(my $jj=1; $jj <= @{$locisref}; $jj++) {
#	    push(@curLi, [ getAlnLociScore($refexons->[$ii -1], $locisref->[$jj -1], $mismV, $penaltyFl, $penaltyFa, $maxRelSc) ] );
	    push(@curLi, [ &{${myScoreFunction}}($refexons->[$ii -1], $locisref->[$jj -1], $mismV, $penaltyFl, $penaltyFa, $maxRelSc, $refscores, $blqueryStI, $bltargetStI) ] );
	    
	}
	push(@scmatrix, \@curLi);
    }

    #print "DESCRL=" . $descrL . "\n";
    my ($imgDPcurX, $imgDPcurY) = ($locSepSize , $locSepSize);
    my $DPmethod = 'filledRectangle';
    
    my $imgPcurX=$locSepSize;
    
    my $mapAreaStr="<div><map name=\"regionID" . $regID . "dpmap\">\n";
    my $mapAreaStrP="<div><map name=\"regionID" . $regID . "dpmapP\">\n";

	## write column numbers / positions for printerDotplot 
	 if($doDPFl) {

		my $qryname=$regDescr;
		my $tmpf=3;
		$tmpf = 1 if(! defined($myfont[- $tmpf]));
			
		$imgDPcurX = $locSepSize;
		$imgDPcurX += $descrL;
		$imgDP->$DPmethod($locSepSize, $locSepSize, $descrL, $eehei, $lblack);
		
		if($printcFl) {      
			$imgPcurX = $locSepSize;
			$imgPcurX += $descrL;
		     #$imgDPp->$DPmethod($locSepSize, $locSepSize, $descrL + scalar(@{$locisref})*($locSepSize + $locSize), 2* $locSepSize+$locSize, $pwhite);
			$imgDPp->$DPmethod($locSepSize, $locSepSize, $descrL, $ehei, $pwhite);
		
			$imgDPp->string($myfont[-$tmpf], int( ($descrL - (length($qryname)*$myfontb[-$tmpf]) - $locSepSize)) +2, $locSepSize , $qryname, $pblack);	
			$imgDPp->string($myfont[-$tmpf], 2 + $locSepSize, $ehei - $locSepSize - $myfonth[-$tmpf] , substr($refDescr, SynBlast::MyUtils::max(0, length($refDescr) - int(($descrL -$locSepSize)/ $myfontb[-$tmpf]) )) , $pblack);	
		}
		
		$imgDP->string($myfont[-$tmpf], int( ($descrL - (length($qryname)*$myfontb[-$tmpf]) - $locSepSize)) +2, $locSepSize , $qryname, $red);	
		$imgDP->string($myfont[-$tmpf], 2 + $locSepSize, $eehei - $locSepSize - $myfonth[-$tmpf] , substr($refDescr, SynBlast::MyUtils::max(0, length($refDescr) - int(($descrL -$locSepSize)/ $myfontb[-$tmpf]) )) , $red);	
					
		for(my $jj=1; $jj <= @{$locisref}; ++$jj) {
			
			my($k1from,$k1to)=split(/_/, [keys %{$locisref->[$jj-1]}]->[0]);
			my $mbasePos=int( ((( ($k1to - $k1from) / 2.) + $k1from) +5000) / 10000.)/100.;
		
			$imgDPp->$DPmethod($imgPcurX, $locSepSize, $imgPcurX+ $locSize, $ehei, $pwhite) if($printcFl);
			$imgDP->$DPmethod($imgDPcurX, $locSepSize, $imgDPcurX+ $locSize, $eehei, $lblack);
				
			if($printcFl) {      
				$imgDPp->stringUp($myfont[-1], $imgPcurX, $ehei - 2, "${mbasePos}", $pblack);
				$imgPcurX += $locSize + $locSepSize;
				$imgPcurX += $colDistSepSizes[$jj-1] if(defined($maxDist) && ($jj<= @colDistSepSizes));
			}				
				
			if(defined($maxDist) && ($jj > 1) && ($colDistances[$jj-2] > $maxDist)) { ## highlight too big distance between columns...
				my $halfSep=int($locSepSize/2);
				$imgDP->$DPmethod($imgDPcurX - $halfSep, $locSepSize, $imgDPcurX + $locSize, $eehei + $locSepSize, $yellow);
				$imgDP->$DPmethod($imgDPcurX, $locSepSize, $imgDPcurX + $locSize, $eehei + $locSepSize, $lblack);
			}
			$imgDP->stringUp($myfont[-1], $imgDPcurX, $eehei - 2, "${mbasePos}", $red);
			$imgDPcurX += $locSize + $locSepSize;
				
				
		}
	}
	
	

    # fill matrix...
    for(my $ii=1; $ii <= @{$refexons}; ++$ii) {
	
	if($doDPFl) {
		$imgDPcurX = $locSepSize;  ##		+ $descrL;
		$imgPcurX = $locSepSize;
	  if(defined($refLidxs) && (index("_" . join("_", @{$refLidxs}) . "_", "_${ii}_") >= 0)) {
	    ### highlight reference locus (focal gene of interest)
	    my $halfSep=int($locSepSize/2);	
	    $imgDP->$DPmethod($imgDPcurX - $halfSep, $imgDPcurY +$eehei - $halfSep, $descrL + $imgDPcurX +  scalar(@{$locisref})*($locSepSize + $locSize) -$locSepSize + $halfSep, $imgDPcurY +$eehei + $locSize + $halfSep, $lgray);
	    $imgDP->$DPmethod($imgDPcurX, $imgDPcurY+$eehei , $descrL + $imgDPcurX +  scalar(@{$locisref})*($locSepSize + $locSize) -$locSepSize, $imgDPcurY +$eehei + $locSize, $black);

	if($printcFl) {
	#    $imgDPp->$DPmethod($imgPcurX- $halfSep, $imgDPcurY +$ehei - $halfSep, $descrL + $imgPcurX +  scalar(@{$locisref})*($locSepSize + $locSize) -$locSepSize + $halfSep, $imgDPcurY +$ehei + $locSize + $halfSep, $pblack);
 	 #   $imgDPp->$DPmethod($imgPcurX, $imgDPcurY + $ehei, $descrL + $imgPcurX+  scalar(@{$locisref})*($locSepSize + $locSize) -$locSepSize, $imgDPcurY +$ehei + $locSize, $plgray);
	  $imgDPp->$DPmethod($imgPcurX- $halfSep, $imgDPcurY +$ehei - $halfSep, $imgDPp->width - 1 - $halfSep, $imgDPcurY +$ehei + $locSize + $halfSep, $pblack);
 	  $imgDPp->$DPmethod($imgPcurX, $imgDPcurY + $ehei, $imgDPp->width - 1  - $locSepSize, $imgDPcurY +$ehei + $locSize, $plgray);
	 }	    
        }
	
	  $imgDPcurX += $descrL;
	$imgPcurX += $descrL;
	
	  ### write query gene name in first "column":
	  if($descrL) {  #== -1) {
	     my $rowDescr="";
	     my $dumy=0;
	     my @qryname=(); #("gene1", "gene2", "gene3");
	     my @crefks = keys %{$refexons->[$ii-1]};
	    (@crefks == 1) or die "Error - There is more than one Hash-Entry for a combined_P_locus ()!\n";
	       my $crefkHash=$refexons->[$ii-1]->{$crefks[0]};
		foreach my $crefkE (keys %{$crefkHash}) {
			foreach my $crefkPE (@{$crefkHash->{$crefkE}}) {
				my $protName=$crefkPE->{'geneName'};
				push(@qryname, $protName);
				$rowDescr .= "; " if($dumy);
				$rowDescr .= $crefkPE->{'protID'};
				$rowDescr .=" [" . $protName . "]" if(length($protName)> 1);
				++$dumy;
				}
		}
         #DataBrowser:browse($refexons->[$ii-1]);
	      my $qryname=join("/", @qryname);
	      $qryname= substr($qryname, 0, $descrFl-1) . ">" if(length($qryname) > $descrFl);
	      $imgDP->$DPmethod($locSepSize, $imgDPcurY+ $eehei, $descrL, $imgDPcurY + $eehei + $locSize, $lblack);
	      $imgDP->string($myfont[-1], int( ($descrL - (length($qryname)*$myfontb[-1]) - $locSepSize)) +2, $imgDPcurY + $eehei , $qryname, $red) if($qryname);	
	      
	      if($printcFl) {
		      $imgDPp->$DPmethod($locSepSize, $imgDPcurY + $ehei, $descrL, $imgDPcurY +$ehei + $locSize, $pwhite);
		      $imgDPp->string($myfont[-1], int( ($descrL - (length($qryname)*$myfontb[-1]) - $locSepSize)) +2, $imgDPcurY + $ehei , $qryname, $pblack) if($qryname);	
	      }
	      
	     $mapAreaStr .= "<area shape=\"rect\" coords=\"" 
		. join(",", ($locSepSize, $imgDPcurY+ $eehei, $descrL, $imgDPcurY + $locSize + $eehei)) 
		. "\" href=\"#regionID_" . $regID . "_p" . $ii . "\" alt=\"p" . $ii . "\"";
		$mapAreaStr .= " title=\"[p${ii}: ${rowDescr}]";
		$mapAreaStr .= "\">\n";
	      
	        $mapAreaStrP .= "<area shape=\"rect\" coords=\"" 
		. join(",", ($locSepSize, $imgDPcurY + $ehei, $descrL, $imgDPcurY + $locSize + $ehei)) 
		. "\" href=\"#regionID_" . $regID . "_p" . $ii . "\" alt=\"p" . $ii . "\""
		. " title=\"[p${ii}: ${rowDescr}]"
		. "\">\n" if($printcFl);
	  }
	  ####
	  
	}
	
	for(my $jj=1; $jj <= @{$locisref}; ++$jj) {
	    my $absSc=$scmatrix[$ii -1]->[$jj -1]->[0];
	    my $relSc=$scmatrix[$ii -1]->[$jj -1]->[2];
	    my $cori = $scmatrix[$ii -1]->[$jj -1]->[3];
	    my $cfrom = SynBlast::MyUtils::formatMyNumber($scmatrix[$ii -1]->[$jj -1]->[4], "'", 3);
	    my $cto = SynBlast::MyUtils::formatMyNumber($scmatrix[$ii -1]->[$jj -1]->[5], "'", 3);
	    
	    my $matchSim = $absSc;
	    
	    if($relAlScFl) {
	    	$matchSim = $relSc if(defined($relSc)); #take relative score as matchSimilarity, if defined (maybe 0!?)
	    }
	
	
	if($doDPFl) {    

	  if(defined($maxDist) && ($jj > 1) && ($colDistances[$jj-2] > $maxDist)) { ## highlight too big distance between columns...
	    my $halfSep=int($locSepSize/2);
	    
	    $imgDP->$DPmethod($imgDPcurX - $halfSep, $imgDPcurY+ $eehei, $imgDPcurX + $locSize, $imgDPcurY + $eehei + $locSize, $yellow);
    	if($printcFl) {
	    my $thickness=$halfSep;
	    $thickness +=	$colDistSepSizes[$jj-2]   if(defined($maxDist) && @colDistSepSizes && ($jj -1 <= @colDistSepSizes));
	    my $bottom=$locSepSize;
	    $bottom=0 if($ii == @{$refexons});
	    $imgDPp->$DPmethod($imgPcurX -$thickness, $imgDPcurY + $ehei, $imgPcurX, $imgDPcurY + $ehei + $locSize + $bottom, $pblack);	    
	}
	
	  }



	  if($matchSim > $mismV) { # match here?
		if(defined($dpColors) && defined($dpColors->[0])) {
		      my $mm= ($maxRelSc / (scalar(@{$dpColors->[0]}) -1));
		     if(defined($relSc)) {
			    my $colIdx = int(($relSc -1)/$mm);
		     	$colIdx = scalar(@{$dpColors->[0]}) -1 if($relSc > $maxRelSc);
			    die "ERROR: ColorIndex out of Range! (alignLocis, ${colIdx}, ${relSc}, ${maxRelSc})\n" unless(($colIdx< scalar(@{$dpColors->[0]})) && ($colIdx>=0));
			    $dred = $imgDP->colorAllocate(@{$dpColors->[0]->[$colIdx]});
  			if(defined($dpColors->[1]) && $printcFl) {
				    $mm= ($maxRelSc / (@{$dpColors->[1]} -1));
				    $colIdx = int(($relSc -1)/$mm);
				    $colIdx = @{$dpColors->[1]} -1 if($relSc > $maxRelSc);
				    $pdred = $imgDPp->colorAllocate(@{$dpColors->[1]->[$colIdx]});
				} 
		      } 
		    } 
		    $imgDP->$DPmethod($imgDPcurX, $imgDPcurY + $eehei, $imgDPcurX + $locSize, $imgDPcurY + $eehei + $locSize, $dred);
        	if($printcFl) {
			$imgDPp->$DPmethod($imgPcurX, $imgDPcurY + $ehei, $imgPcurX + $locSize, $imgDPcurY +$ehei + $locSize, $pdred);
			$imgDPp->$DPmethod($imgPcurX-1, $imgDPcurY + $ehei -1, $imgPcurX + $locSize +1, $imgDPcurY +$ehei + $locSize +1, $pdred) 
				if( defined($relSc) && ($relSc > $maxRelSc));
		}


	} else { ## no match
		    $imgDP->$DPmethod($imgDPcurX, $imgDPcurY + $eehei, $imgDPcurX + $locSize, $imgDPcurY + $eehei + $locSize, $gray);

	if($printcFl) {
		    $imgDPp->$DPmethod($imgPcurX, $imgDPcurY + $ehei, $imgPcurX + $locSize, $imgDPcurY + $ehei + $locSize, $pwhite);
	    }
	    
}
	    
	   	if($cori) { ### at least one misorientated query-hit  in current p-locus?  draw a shadow...
	   		$imgDP->$DPmethod($imgDPcurX, $imgDPcurY + $eehei, $imgDPcurX + int(1/8. * $locSize), $imgDPcurY + $eehei + $locSize, $lblack);
			$imgDP->$DPmethod($imgDPcurX, $imgDPcurY+ $eehei, $imgDPcurX + $locSize, $imgDPcurY + $eehei + int(1/8. * $locSize), $lblack);
	        	if($printcFl) {
				$imgDPp->$DPmethod($imgPcurX, $imgDPcurY + $ehei, $imgPcurX + int(1/8. * $locSize) -1, $imgDPcurY + $ehei + $locSize, $pwhite);
				$imgDPp->$DPmethod($imgPcurX, $imgDPcurY + $ehei, $imgPcurX + $locSize, $imgDPcurY +$ehei + int(1/8. * $locSize) -1, $pwhite);
	        	}
	   	} 

	    $mapAreaStr .= "<area shape=\"rect\" coords=\"" 
		. join(",", ($imgDPcurX, $imgDPcurY + $eehei, $imgDPcurX + $locSize,  $imgDPcurY + $eehei + $locSize)) 
		. "\" href=\"#regionID_" . $regID . "_p" . $ii . "\" alt=\"p" . $ii . ", q" . $jj . "\"";
		$mapAreaStr .= " title=\"[p${ii}-q${jj}]";
		
		$mapAreaStrP .= "<area shape=\"rect\" coords=\"" 
		. join(",", ($imgPcurX, $imgDPcurY + $ehei, $imgPcurX + $locSize,  $imgDPcurY + $locSize + $ehei)) 
		. "\" href=\"#regionID_" . $regID . "_p" . $ii . "\" alt=\"p" . $ii . ", q" . $jj . "\""
		. " title=\"[p${ii}-q${jj}]" if($printcFl);
		
		if($relAlScFl && defined($relSc)) {
			$mapAreaStr .= ": rscore=${relSc}, score=${absSc}, coord=${cfrom}-${cto}" if(${absSc}> $mismV);  #, ori=${cori}";
			$mapAreaStrP .= ": rscore=${relSc}, score=${absSc}, coord=${cfrom}-${cto}" if($printcFl && ($absSc> $mismV));  #, ori=${cori}";
		} else {
				$mapAreaStr .= ": score=${absSc}" 
				            . (defined($relSc) ? ", rscore=${relSc}" : "") 
				            . ", coord=${cfrom}-${cto}" if(${absSc}> $mismV);  #, ori=${cori}";
				$mapAreaStrP .= ": score=${absSc}" 
				            . (defined($relSc) ? ", rscore=${relSc}" : "") 
				            . ", coord=${cfrom}-${cto}" if($printcFl && ($absSc> $mismV));  #, ori=${cori}";
			} 
		$mapAreaStr .= "\">\n";
		$mapAreaStrP .= "\">\n" if($printcFl);


	    $imgDPcurX+=$locSize + $locSepSize;
	    $imgPcurX+=$locSize + $locSepSize;
	    $imgPcurX+= $colDistSepSizes[$jj-1] if(defined($maxDist) && @colDistSepSizes && ($jj <= @colDistSepSizes));
	    
	  } ## ifDPflag
	  
	  
	    my $top = $matrix[$ii -1]->[$jj] + $gapsimQ; # gap-insert in target-loci-sequence # should be free or with lower cost (if some genes of informant sequence (ploci) are not found, this can be ignored
	    my $left = $matrix[$ii]->[$jj -1] + $gapsimP; # gap-insert in query-exon-sequence # cost should be larger because of possible dissimilarity of order of genes
	    my $diag = $matrix[$ii -1]->[$jj -1] + $matchSim;
	    push(@{$matrix[$ii]}, SynBlast::MyUtils::max($diag, $top, $left));
	    #push(@{$matrix[$ii]}, max( 0, max( max($diag, $top), $left)));
	    # $matrix[$ii]->[$jj] = max( max($diag, $top), $left);
	}
	$imgDPcurY += ($locSize + $locSepSize) if($doDPFl);
    }
    
    $mapAreaStr .= "</map></div>\n";
    $mapAreaStrP .= "</map></div>\n" if($printcFl);
    
    
    ## do backtracking...
    my @alnStr= ();
    #Score of Needleman-Wunsch-Global-Aln is Maximum of matrix-borders (as endgaps should be free (no costs) here...)
    my ($lastI, $lastJ, $alnScore) = getMatrixSidePos(\&SynBlast::MyUtils::max, \@matrix);

    print "Matrix-Maximum data:\n" . join("\t", ( $lastI, $lastJ, $alnScore)) . "\n" if($verboseFl);

    # do backtracking (endfreeGaps) 
    my $bstartI=$lastI; #@{$refexons};
    my $bstartJ=$lastJ; #@{$locisref};
    # do backtracking, get one optimal alignment-String for calculated Score... 
    if($lastI == @{$refexons}) {
	## ploci are already complete, just add the qloci aligned to gaps (gape)
	for(my $rJ=$lastJ+1; $rJ <= @{$locisref}; $rJ++) {
	    push(@alnStr, [ ( "---", "q" . $rJ, "gape") ] );
	}
    } else {
	## qloci are already complete, just add the ploci aligned to gaps (gape)
	for(my $rI=$lastI+1; $rI <= @{$refexons}; $rI++) {
	    push(@alnStr, [ ( "p" . $rI, "---", "gape") ] );
	}
    }
    # get rest of alignmentStr starting at maximum score position in a matrix side
    
    my @alnIdxA=();
    traceback($bstartI, $bstartJ, \@alnStr, \@matrix, $gapsimP, $gapsimQ, $mismV, \@scmatrix, "---", \@alnIdxA, $relAlScFl);

    #put out aln-String
    SynBlast::MyUtils::printListList(\@alnStr, "\t", "\n") if($verboseFl);
    
    
    my @DPdata=();
    my @DPdatap=();

 if($doDPFl) {
	 my $cumuStart=0;
	 my @cumuSepSizes=();
	 for(my $dsli=0; $dsli< @colDistances; ++$dsli) {
		$cumuStart+=$colDistSepSizes[$dsli] if(@colDistSepSizes);
		push(@cumuSepSizes, $cumuStart);
	}
	
    foreach my $cDot (@alnIdxA) {
	my ($ppp, $qqq, $sc, $relSc) = @{$cDot};
	my $imgDPcurX = $descrL + $locSepSize + ($locSize + $locSepSize)*($qqq -1) + int($locSize/3 +0.5);
	
	my $imgPcurX = $imgDPcurX;
	$imgPcurX+= $cumuSepSizes[$qqq-2] if($qqq>1);
	
	my $imgDPcurY = $locSepSize + ($locSize + $locSepSize)*($ppp -1) + int($locSize/3 +0.5);
	$imgDP->$DPmethod($imgDPcurX, $imgDPcurY + $eehei, $imgDPcurX + int($locSize/3), $imgDPcurY + $eehei + int($locSize/3), $black);
	if($printcFl) {
		$imgDPp->$DPmethod($imgPcurX, $imgDPcurY + $ehei, $imgPcurX + int($locSize/3), $imgDPcurY + $ehei + int($locSize/3), $pwhite);
	}
	
    }

    my $DPpngStr= $imgDP->png;
    @DPdata = (\$DPpngStr, \$mapAreaStr);

    if($printcFl) {
	my $DPpngStrp=$imgDPp->png;
	@DPdatap = (\$DPpngStrp, \$mapAreaStrP) if($DPpngStrp);
    }
    
#    my $DPpngStr= $imgDP->jpeg(50);
}
  
  return($alnScore,\@alnStr, \@DPdata, \@alnIdxA, \@DPdatap);
    
    



    sub getMatrixPos {
	my ($methodRef, $listref) = @_;
	my ($idxEins, $idxZwei, $ovalue) = ( undef, undef, undef);
  	for(my $lEins=0; $lEins< @{$listref}; $lEins++) {
 	    for(my $lZwei=0; $lZwei< @{$listref->[$lEins]}; $lZwei++) {
		my $cvalue=$listref->[$lEins]->[$lZwei];
		($idxEins, $idxZwei, $ovalue) = ($lEins, $lZwei, $cvalue) if ($cvalue == &{$methodRef}($cvalue, $ovalue));
	    }
	}
#	&{$methodRef}(,)
	return($idxEins, $idxZwei, $ovalue);
    }
    
    
    sub getMatrixSidePos {
	my ($methodRef, $listref) = @_;
	my ($idxEins, $idxZwei, $ovalue) = ( undef, undef, undef);
  	for(my $lEins=0; $lEins< @{$listref} -1; $lEins++) {
 	    my $cvalue=$listref->[$lEins]->[@{$listref->[$lEins]} -1];
	    ($idxEins, $idxZwei, $ovalue) = ($lEins, (@{$listref->[$lEins]} -1), $cvalue) if ($cvalue == &{$methodRef}($cvalue, $ovalue));
	}
	for(my $lZwei=0; $lZwei< @{$listref->[@{$listref} -1]}; $lZwei++) {
	    my $cvalue=$listref->[@{$listref} -1]->[$lZwei];
	    ($idxEins, $idxZwei, $ovalue) = (@{$listref} -1, $lZwei, $cvalue) if ($cvalue == &{$methodRef}($cvalue, $ovalue));
	}
	
#	&{$methodRef}(,)
	return($idxEins, $idxZwei, $ovalue);
    }
    
    
    sub getAlnLociScore {
	my ($curPloci, $curLoci, $mismV, $penaltyFl, $penaltyFa, $maxRelSc, $refscores) = @_;
	## returns (score, \successStr, relativeScore)
	## successStr is reference on string showing the number of correct matches (incl. orientation, if flag set)
	$penaltyFl=0 if(! defined($penaltyFl));
	
	# $mismV=-1e2 if(! defined($mismV));
	# $penaltyFa=10 if(! defined($penaltyFa));  ###(10 percent of match score is used in case of contrary orientation)
	
	my @lkey= keys %{$curLoci};
	(@lkey == 1) or die "Error: There are more than one Hash-Entries for a putative loci (getAlnLociScore)\n";
	my @plkey= keys %{$curPloci};
	(@plkey == 1) or die "Error: There are more than one Hash-Entries for a query loci (ploci) (getAlnLociScore)\n";
	
	####
	### matchScore is the Sum of all matching queries out of current ploci in current qloci
	### relMatchScore is matchScore divided by Sum of all matching queries's selfScore
	my $matchV=0;
	my $rmatchV=0;
	my $cTotal=0;
	my $cSuccess=0;
	
	my $minPos=undef;  ## for displaying the interval of matching part of the qLocus
	my $maxPos=undef;
	
	
	my $cQsuccess=0;
	foreach my $curRefk (keys %{$curPloci->{$plkey[0]}}) {
	    
	    my $maxCEV=0;
	    my $maxCErV=0;
	    my $mCEVtotal=0;
	    my $mCEVsucc=0;
	    foreach my $cE (@{$curPloci->{$plkey[0]}->{$curRefk}}) {
		# go through all members of the current loci and keep only the best result
		# (as overlapping/matching ploci are assumed to be the same)
		my $cE_matchV=0;
		my $cE_rmatchV=0;
		my $cE_succ=0;
		my $cE_total=0;
		

		my $curRefID = $cE->{'protID'};
		
		## old version: selfScore entry from infoFile is used for relativeScore calculation
		##my $curRefsSc = (exists($cE->{'selfScore'}) ? int($cE->{'selfScore'}) : undef);
		
		### now: the globally-best hit is used as selfScore (as determined within the same program run)
		my $curRefsSc=undef;
		$curRefsSc= $refscores->{$curRefID}->[0] if(exists($refscores->{$curRefID}) ); ## && (@{$refscores->{$curRefID}}));					 
					 
					 
		

		if(exists($curLoci->{$lkey[0]}->{$curRefID})) {
		    #match?

		    foreach my $cL (@{$curLoci->{$lkey[0]}->{$curRefID}}) {
			++${cE_total};
			if(($penaltyFl) && misOrientation($cE, $cL, $penaltyFl)) {
			    $cE_matchV += (($cL->{'score'}) * $penaltyFa / 100.0);
			    $cE_rmatchV += ((($cL->{'score'}) * $penaltyFa / 100.0) / ${curRefsSc} * $maxRelSc) if(${curRefsSc});
			} else {
			    $cE_matchV += ($cL->{'score'}); 
			    $cE_rmatchV += (($cL->{'score'})/${curRefsSc} * $maxRelSc) if(${curRefsSc});
			    ++${cE_succ};
			}
			$minPos=SynBlast::MyUtils::min($minPos, $cL->{'startPos'});
			$maxPos=SynBlast::MyUtils::max($maxPos, $cL->{'endPos'});
    }
		  
		}
		
		
		# check if current matchV is best
		if(${cE_matchV} > $maxCEV) {
		    $maxCEV= ${cE_matchV};
		    $maxCErV= ${cE_rmatchV};
		    $mCEVtotal=${cE_total};
		    $mCEVsucc=${cE_succ};
		}
	    }
	    $matchV+=$maxCEV;
	    $rmatchV+=$maxCErV;
	    $cTotal+=$mCEVtotal;
	    $cSuccess+=$mCEVsucc;
	    ++$cQsuccess if($maxCEV);
	}
	
	if($cQsuccess > 1) {
	    $matchV = $matchV / $cQsuccess;
	    $rmatchV = $rmatchV / $cQsuccess;
	}

	
	$matchV=int($matchV); # because of errors in backtracking otherwise
	$rmatchV=int($rmatchV); # because of errors in backtracking otherwise
	
	if($cTotal) {
	    my $sStri= ("+" x $cSuccess);
	    $sStri= $sStri . "/" . ("+" x $cTotal) if($penaltyFl);
	    return($matchV, \$sStri, $rmatchV, ($cTotal - $cSuccess), $minPos, $maxPos) if($matchV >0);
	    ## (cTotal - cSuccess) is the number of misorientated hits in current p-q-assignment
	}
	return($mismV, undef, undef, undef, undef, undef);
    }
    
    sub misOrientation {
	my ($plociEntry, $qlociEntry, $penFl) = @_;
	return(0) if(! $penFl);
	
	my $plociO=$plociEntry->{'ori'}; ## -1 or 1 for query-orientation on slice (strand_field)
	#(see infofile; fields are: cdsStart_cdsEnd_ProteinID_TranscriptID_GeneID_Strand_PeptidLength_MeanCdsPos_ExternalName
	
	my $qlociO=$qlociEntry->{'ori'}; 
	#fields are: start_end_BlastHitsNumber_orient_meanBlastHitPos_score
	# qlociOri is 1 if protStart < protEnd (the same as on DNA_target_side)
	# 
	if(($penFl == 2)) {
	    # qlocis in decreasing order...(revCombLocis is used) -> revert also the qlociOrientation
	    $qlociO= ($qlociO > 0) ? -1 : 1;
	}
	
	#if(($plociO > 0)) {
	#    # ploci-Protein orientation on slice is same as slice (ploci-sequence)
	#    return(0) if(
	#}
	return($plociO != $qlociO);
    }
    
    

    sub traceback {
	my ($curI, 
	    $curJ, 
	    $alnStrA, 
	    $matRef, 
	    $gapsimP,
	    $gapsimQ,
	    $mismV,
	    $scmatRef,
	    $gapchar, 
	    $alnIdxA,
	    $relAlScFl
	      ) = @_;
	my $gapType="gap";
	
	while($curI || $curJ) { # as long as neither curI nor curJ is zero...
	    my $matchS = 0;
			my $relMatchS = 0;
			my $absMatchS = 0;
	    if(($curI >0) && ($curJ >0) ) {
		$absMatchS = $scmatRef->[$curI -1]->[$curJ -1]->[0];
		$relMatchS = $scmatRef->[$curI -1]->[$curJ -1]->[2];
		$matchS = $absMatchS;
		if($relAlScFl) {
			$matchS = $relMatchS if(defined($relMatchS));
		}
		
		if ( ($matRef->[$curI]->[$curJ] - $matRef->[$curI -1]->[$curJ -1]) == $matchS) {
		    # diagonal way possible?
		    my $matchChar="-";
		    $matchChar= ${$scmatRef->[$curI -1]->[$curJ -1]->[1]} if( $matchS > $mismV); # match?  
		    splice(@{$alnStrA}, 0, 0, [ ( "p" . $curI, "q" . $curJ, $matchChar) ] );
		    push(@{$alnIdxA}, [ ($curI, $curJ, $absMatchS, $relMatchS) ]);
		                  # stores matching loci indices and score and relScore
		    --$curI;
		    --$curJ;
		    next;
		}
	    }

	    if( ($curI > 0) && ( ($matRef->[$curI]->[$curJ] - $matRef->[$curI -1]->[$curJ]) == $gapsimQ) ) {
		# gap insert in target-loci-sequence possible...
		splice(@{$alnStrA}, 0, 0, [( "p" . $curI, $gapchar, $gapType . "q") ] );
		--$curI;
		next;
	    }
	    
	    if( ($curJ > 0) && ( ($matRef->[$curI]->[$curJ] - $matRef->[$curI]->[$curJ -1]) == $gapsimP ) ) {
		# gap insert in query-sequence possible...
		splice(@{$alnStrA}, 0, 0, [( $gapchar, "q" . $curJ, $gapType . "p") ] );
		--$curJ;
		next;
	    }
	    
	    
	    if( ($curJ > 0) && ( $matRef->[$curI]->[$curJ -1] == 0 ) ) {
		# startgap insert in query-sequence possible...
		splice(@{$alnStrA}, 0, 0, [( $gapchar, "q" . $curJ, "gaps") ] );
		--$curJ;
		next;
	    }
	
	    if( ($curI > 0) && ( $matRef->[$curI -1]->[$curJ] == 0) ) {
		# startgap insert in target-loci-sequence possible...
		splice(@{$alnStrA}, 0, 0, [( "p" . $curI, $gapchar, "gaps") ] );
		--$curI;
		next;
	    }

	    die "FATAL Error: traceback!\n" 
		. $curI . ", " . $curJ . "\n"
		. ($matRef->[$curI]->[$curJ] - $matRef->[$curI -1]->[$curJ -1]) . "\n"
		. $matchS . "\n"
		. ($matRef->[$curI]->[$curJ] - $matRef->[$curI -1]->[$curJ]) . "\n"
		. ($matRef->[$curI]->[$curJ] - $matRef->[$curI]->[$curJ -1]) . "\n"
		. $gapsimP 
		. $gapsimQ
		. "\n";
	}
	
	return(1);
    }
 

    
    
}  ## end  of sub alignLocis







##########################################################################
# combines target loci into non-overlapping q-loci, by grouping overlapping chained-HSP intervals (of different queryID)
# according to start/end or mean HSP positions
##########################################################################
sub combineLoci2List {
    my ($locilist, $outref, $sFlag, $qLflag) = @_;
    ### if qLflag is set, locilist is assumed to be in q-loci-format, and 
    ### overlapping entries for a range with the same protQueryID are reduced to the best-scoring one

    @{$outref}=();
    $sFlag=0 if(! defined($sFlag));

    $qLflag=0 if(! defined($qLflag));
    

    
    for(my $li=@$locilist -1;$li>=0; $li--) {
	my $curRangeEnd=0;
	my $curRangeStart=1e11;
	my %locis=(); # to be a hash{key:DNA-range of combined loci} with value=hash{key:QueryID} with value=query-lociRanges_blastHitsNumber_orientation_lociScore


	for my $ckey (keys %{$locilist->[$li]}) {
	    my $cvalues=$locilist->[$li]->{$ckey}; # is now a hashref
	    if(! $sFlag) {
		$curRangeEnd=SynBlast::MyUtils::max($curRangeEnd, $cvalues->{'endPos'});
		$curRangeStart=SynBlast::MyUtils::min($curRangeStart, $cvalues->{'startPos'});
	    } else {
		$curRangeEnd=SynBlast::MyUtils::max($curRangeEnd, $cvalues->{'meanPos'});
		$curRangeStart=SynBlast::MyUtils::min($curRangeStart, $cvalues->{'meanPos'});
	    }
	}
	my $rli= $li -1;
	my $rRangeEnd=$curRangeStart;
	my @idxToAdd=();
	while(($rli >= 0) && ($rRangeEnd >= $curRangeStart)) { # as long as last block did overlap
	    my $rRangeEnd=0;
	    my $rRangeStart=1e11;
	    # obtain range for block..
	    for my $ckey (keys %{$locilist->[$rli]}) {
		my $cvalues=$locilist->[$rli]->{$ckey};
		if(! $sFlag) {
		    $rRangeEnd=SynBlast::MyUtils::max($rRangeEnd, $cvalues->{'endPos'});
		    $rRangeStart=SynBlast::MyUtils::min($rRangeStart, $cvalues->{'startPos'});
		} else {
		    $rRangeEnd=SynBlast::MyUtils::max($rRangeEnd, $cvalues->{'meanPos'});
		    $rRangeStart=SynBlast::MyUtils::min($rRangeStart, $cvalues->{'meanPos'});
		}
	    }
	    ($rRangeEnd <= $curRangeEnd) or die "Error: loci-list not sorted correctly! (combineLoci2List)\n";
	    
	    if(SynBlast::MyUtils::getOverlapSize($curRangeStart, $curRangeEnd, $rRangeStart, $rRangeEnd) > 0) {
		# overlaping block found -> add
		push(@idxToAdd, $rli);
		$curRangeStart=SynBlast::MyUtils::min($curRangeStart, $rRangeStart);
	    }
	    $rli--;
	}


	my %curRangeHash=();
	# put in the original block-members
	## only one entry per protQueryID (q-loci-case) and range should already be the result of calcSplAln
	##
	for my $ckey (keys %{$locilist->[$li]}) {
	    push( @{$curRangeHash{$ckey}}, { %{$locilist->[$li]->{$ckey}} } );
	}
	
	# add the members of overlapping blocks...
	


	## if qLflag set, only the best scoring hit per range and protQueryID is kept... !
	if($qLflag) {
	    my $tdumy=$li;
	    my $needToRecalc=0;
	    for(my $ji=0; $ji< @idxToAdd; $ji++) {
		#test on error?
		$tdumy--;
		($tdumy == $idxToAdd[$ji]) or die "Error: idxToAdd contains invalid entries (combineLoci2List)\n";
		# add overlapping block to curRangeHash
		for my $ckey (keys %{$locilist->[$idxToAdd[$ji]]}) {
		    if(exists($curRangeHash{$ckey})) {
			## already an entry for the same protID-range-combination???
			#if($curRangeHash{$ckey}->[0]->[-1] < $locilist->[$idxToAdd[$ji]]->{$ckey}->[-1]) {
			if($curRangeHash{$ckey}->[0]->{'score'} < $locilist->[$idxToAdd[$ji]]->{$ckey}->{'score'}) {
			    ## new entry has better score for this combination?
			    # replace
			    $needToRecalc++;
			    @{$curRangeHash{$ckey}} = ( { %{$locilist->[$idxToAdd[$ji]]->{$ckey}} } );
			} 
		    } else {
			## no entry yet for this protID-range-combination???
			## push as first entry
			push( @{$curRangeHash{$ckey}},  { %{$locilist->[$idxToAdd[$ji]]->{$ckey}} } );
		    }
		}
	    }
	    
	    if($needToRecalc) {
		## calculate the range boundaries again, as there might have been a change
		$curRangeStart=1e11;
		$curRangeEnd=0;
		for my $ckey (keys %curRangeHash) {
		    $curRangeStart=SynBlast::MyUtils::min($curRangeStart, $curRangeHash{$ckey}->[0]->{'startPos'});
		    $curRangeEnd=SynBlast::MyUtils::max($curRangeEnd, $curRangeHash{$ckey}->[0]->{'endPos'});
		}
		print STDERR "\n >recalculated the boundaries for current block, as changed\n";
	    }
	}
	
	else {  ### several Hits for the same Range and ckey allowed (ploci-case, to support double names etc)
	    my $tdumy=$li;
	    for(my $ji=0; $ji< @idxToAdd; $ji++) {
		#test on error?
		$tdumy--;
		($tdumy == $idxToAdd[$ji]) or die "Error: idxToAdd contains invalid entries (combineLoci2List)\n";
		# add overlapping block to curRangeHash
		for my $ckey (keys %{$locilist->[$idxToAdd[$ji]]}) {
		    push( @{$curRangeHash{$ckey}},  { %{$locilist->[$idxToAdd[$ji]]->{$ckey}} } );
		}
	    }
	}
	
	$locis{$curRangeStart . "_" . $curRangeEnd} = \%curRangeHash; # hash-use only for storing the block-range-borders
	$li= $li - @idxToAdd; # skip already included blocks...
	
	push(@{$outref}, \%locis);
    }
} 







##########################################################################
#  sorting target loci according to End/Mean position
##########################################################################
sub qlociSort2List {
    my ($lociScRef, $outref, $sFlag) = @_;
    @{$outref}=();
    
    $sFlag=0 if(! defined($sFlag));

    my %locis=(); # to be later a sortable hash{key:DNA-end-str} with value=hash{key:QueryID} with value=array(lociRanges,lociScore...)
    foreach my $curQ ( keys %{$lociScRef} )  {
	for(my $li=0; $li < @{$lociScRef->{$curQ}}; ++$li) {
	    my $curSortKey=int($lociScRef->{$curQ}->[$li]->{'endPos'}); 
	    # sort by dna-end-pos (normal case, sFlag not set)

	    if($sFlag) {
		$curSortKey=int($lociScRef->{$curQ}->[$li]->{'meanPos'}); # MeanBlHitPos as sort index
	    }
	    
	    if(defined($locis{$curSortKey}->{$curQ})) {
		die "Error: another Loci-Pair for the same Query and (mean/interv)-BlastHit DNA-position already existing! (qlociSort2List)\n";
	    }
	    
	    $locis{$curSortKey}->{$curQ}= { %{$lociScRef->{$curQ}->[$li]} };
	}
    }

    my @sortedKeys= (sort {$a <=> $b} (keys %locis)); 
    for(my $li=0; $li<@sortedKeys; $li++) {
	push(@{$outref}, { %{$locis{$sortedKeys[$li]}} } );
    }
}





sub checkLociLengths {
    my ($lociRef, $exonl, $queryIdxs, $factor, $icolIdx) = @_;
    if (! defined($factor)) { 
	print "loci larger than (reference) are:\n"; $factor=1; 
    } else { print "loci larger than ${factor}x(reference) are:\n"; }
    
    foreach my $curP (keys %{$lociRef}) {
	my $exonRef=  \@{$exonl->[$queryIdxs->{$curP}]};
	my $oriL=abs($exonRef->[$icolIdx->{'endPos'}] - $exonRef->[$icolIdx->{'startPos'}]) +1;
	for(my $lid=0; $lid< @{$lociRef->{$curP}}; $lid++) {
	    my $curL= abs($lociRef->{$curP}->[$lid]->{'endPos'} - $lociRef->{$curP}->[$lid]->{'startPos'}) +1;
	    if ($factor * $oriL < $curL) {
		print $curP . "\t" . $lid . "\t";
		SynBlast::MyUtils::printHash($lociRef->{$curP}->[$lid], "=", ",");
		print "\t" . $curL . "\t(" . $oriL . ")\n";
	    }
	}
    }
}






##########################################################################
# HSP chaining procedure using alignLocis to match consistent orders of interval sets
##########################################################################
sub calcSplAln2 {
 my ($lref, #$querylref, 
     $resultRef, $queryID, $maxSize, $scutoff, $maxOvl, $verboseFl, $queryL, $queryStartI, $targetStartI, $lenI, $preLcds, $posDiffFac) = @_;

$posDiffFac=3 unless(defined($posDiffFac));
 $verboseFl=0 if(! defined($verboseFl));
 $queryStartI=6 if(! defined($queryStartI));
 $targetStartI=8 if(! defined($targetStartI));
 $lenI=3 if(! defined($lenI));

# assumption!!!: lref already sorted ascending acc. to target/DNA-end positions !!!
 my $dummy=0;
 for(my $li=0; $li< @{$lref}; ++$li) {

   if($lref->[$li]->[$targetStartI +1] >= $dummy) {
     $dummy=$lref->[$li]->[$targetStartI +1];
   } else {
     die "ERROR: lref not sorted asc. in target-end-position!\n";
   }
 }
 return(0) if(! @{$lref});
 
 # build target-sets of blasthit-intervals, separated by orientation
 my $mOvl = 0;
 my $targetLi;
 @{$targetLi} = reverse @{$lref};
 
 my $targetLiPos=filterOrientation($targetLi, 1, $queryStartI, $targetStartI);
 my $targetLiNeg=filterOrientation($targetLi, 0, $queryStartI, $targetStartI);
 
 my $targetRegionsPos = buildBlastHitSet($targetLiPos, $targetStartI +1,  $mOvl, $maxSize, 0, $verboseFl, $queryStartI, $targetStartI);
 my $targetRegionsNeg = buildBlastHitSet($targetLiNeg, $targetStartI +1,  $mOvl, $maxSize, 0, $verboseFl, $queryStartI, $targetStartI);
 
 my $targetArr;
 @{$targetArr}= ( $targetRegionsPos, $targetRegionsNeg );


### now perform alignment of every targetRegion to the one queryRegion
 my $myAlnScoreFun = \&getSplAlnScore;
 my $oCountr=0;
 foreach my $targetRegions (@{$targetArr}) {
   foreach my $targR (@{$targetRegions}) {
     ### do alignment twice (for both directions)
     
     # get the corresponding queryRegion
     my $cQueryHash;
     %{$cQueryHash} =( 'unsorted' => [ () ] );
     foreach my $tmpLi (@{$targR}) {
       my @kes = keys %{$tmpLi};
       die "Error: only one key allowed per array here!\n" unless(@kes == 1);
       foreach my $tmpAr (@{$tmpLi->{$kes[0]}}) {
	 push(@{$cQueryHash->{'unsorted'}}, $tmpAr);
       }
     }
     # now sort target specific hits for query positions...
     my $queryHashRef;
     %{$queryHashRef}=();
     if(! (SynBlast::MySyntenyUtils::lists2sortedLists($cQueryHash, $queryHashRef, 1, $queryStartI +1, 0, $queryID, 0, undef, $queryStartI, $targetStartI) > 0)) {
       print STDERR "error: protHits-Hash contains no data (at least for query \"${queryID}\"! (calcSplAln2)\n";
       return(undef);
     }
     my @chrs = keys %{$queryHashRef};
     die "ERROR: Only one type of chromosome is expected for target-specific query-List!\n" unless(@chrs == 1);
     
     my $queryLi;
     @{$queryLi} = reverse @{$queryHashRef->{$chrs[0]}};
     my $queryRegions =  buildBlastHitSet($queryLi, $queryStartI +1, $maxOvl, $maxSize, 0, $verboseFl, $queryStartI, $targetStartI);
     die "ERROR: Only one group is expected for query-List!\n" unless(@{$queryRegions} == 1); 

  ### now perform alignment
     my $tR_queryR;
     
     if($oCountr) { # negative Orientated hits?
       @{$tR_queryR} = reverse @{$queryRegions->[0]};
     } else {
       @{$tR_queryR} = @{$queryRegions->[0]};
     }
     
     my ($alnScore
	 ,$alnStrRef
	 , $alnDPref
	 , $alnIdxR, $alnDPrefp)=alignLocis($targR       #$refexons
				, $tR_queryR #$locisref
				, -1         #$gapsimP
				, -1         #$gapsimQ
				, -1000      # $mismV
				, $verboseFl #
				, 1          # $penaltyFl
				, 1          # $penaltyFa
				, "noRegion" # $regID
				, undef      # $dpColors
				, undef ## printcFl
				, $myAlnScoreFun #
				, 0          # $relAlScFl
				, undef       # $refLidxs
				, 0      # $doDPFl -> do not create png string for dot plot
				, undef ## no maxDistance param
				, $queryStartI 
				, $targetStartI 
				);
     
     
#     DataBrowser::browse($alnStrRef);
     
     #extract actuall hits corresponding to the aln-matches
     my @selEntries=();
     my $finalOri=undef;
     my $scMax=0;
     for(my $li=0; $li< scalar(@{$alnStrRef}); ++$li) { # go through alnStr
       my $cP = $alnStrRef->[$li]->[0]; # e.g. p1 or ---
       my $cQ = $alnStrRef->[$li]->[1]; # e.g. q4 or ---
       my $cplo = undef;
       $cplo = $targR->[int(substr($cP, 1))-1] if($cP ne "---");    #ploci-data
       my $cqlo =undef;
       $cqlo = $tR_queryR->[int(substr($cQ, 1))-1] if($cQ ne "---");  # qloci-data 
       if(defined($cplo) && defined($cqlo)) {
	 #find maximum pairwise match again...
	 my($matchV, $bestHit, $cMatchOri, $cMatch, $cOri)=getMaxPairwiseMatch($cplo, $cqlo, 1, 1,  $queryStartI, $targetStartI);
	 push(@selEntries, $bestHit);
	 $scMax+=$matchV;
	 die "ERROR: NO ori available!\n" unless(defined $cOri);
	 if($matchV == 0) { print "ERROR: MaxPairwiseMatch has value zero!!!\n"; }
	 $finalOri= $cOri if(! defined($finalOri));
	 die "ERROR: Ori changed!\n" unless($finalOri == $cOri);
	 
       }
     }
     
     
     next if($scMax < $scutoff);
          
     # now all hits are defined as one assembled blasthit...
     @selEntries = reverse @selEntries;
     my @curRangeField=getRangeField(\@selEntries, $queryStartI, $targetStartI, $lenI, $preLcds, $posDiffFac);
     
     my %hcontent=( 'score' => $scMax
		    , 'hits' => \@selEntries
		    , 'ori' => $finalOri
		    , 'startPos' => $curRangeField[0]
		    , 'endPos' => $curRangeField[1]
		    , 'hitsCount' => $curRangeField[2]
		    , 'meanPos' => $curRangeField[3],
		    , 'hitsLength' => $curRangeField[4]
		    , 'queryCov' =>  ( int($queryL) ? ("" . (int(( $curRangeField[5] / $queryL * 1000) + 0.5) / 10) . "%") : ("n.a."))
		    , 'chr' => $curRangeField[6]
		    , 'meanIdentity' => ( "" . $curRangeField[7] ."%")
		    , 'CDSannot' => $curRangeField[8]
		    , 'queryID'  => $selEntries[0]->[0]
		    
		   ,  'uncoveredAA1' => ( int($queryL) ?  ( ($curRangeField[9]->[0] -1) . "(" . (int( ( ($curRangeField[9]->[0] -1) / $queryL * 1000.) + 0.5) /10.) . "%)" ) : ( ($curRangeField[9]->[0] -1) ))
		   , 'uncoveredAA2'  => ( int($queryL) ?  ( ($queryL - $curRangeField[9]->[1]) . "(" . (int( ( ($queryL - $curRangeField[9]->[1]) / $queryL * 1000.) + 0.5) /10.) . "%)" ) : ("n.a.")) #( ($queryL - $curRangeField[9]->[1]) . "aa"))
	
		    #   , 'org' => ...
		    );
     
     ## put here check on overlap with prior locis for queryID
     if(checkLociOverlapOK(\@{$resultRef->{$queryID}}, \@curRangeField)) {
       push(@{$resultRef->{$queryID}}, \%hcontent);
     } else {
       print "Putative Loci \"" . join("_", @curRangeField) . "\t${scMax}\" is overlapping with a higher scoring loci -> skipping this entry...\n" if($verboseFl);
     }
   }
   ++$oCountr;
 }
 
 print "END\n\n" if($verboseFl);
 
 
 
 
 sub getSplAlnScore {
   my ($curPloci, $curLoci, $mismV, $penaltyFl, $penaltyFa, $maxRelSc, $refscores, 
		$queryStartI, $targetStartI) = @_;
   ## returns (score, \successStr, relativeScore)
   ## successStr is reference on string showing the number of correct matches (incl. orientation, if flag set)
   $penaltyFl=0 if(! defined($penaltyFl));
   $queryStartI=6 if(! defined($queryStartI));
   $targetStartI=8 if(! defined($targetStartI));
   
   ####
   ### matchScore is the Sum of all matching position intervalls out of current ploci in current qloci
   ### relMatchScore is matchScore divided by Sum of all matching queries's selfScore
   
   my($matchV, $bestHit, $cMatchOri, $cMatch, $cOri)=getMaxPairwiseMatch($curPloci, $curLoci, $penaltyFl, $penaltyFa,  $queryStartI, $targetStartI);
   
   $matchV=int($matchV); # because of errors in backtracking otherwise...?
   if($cMatch) {
     my $sStri= ("+" x $cMatchOri);
     $sStri= $sStri . "/" . ("+" x $cMatch) if($penaltyFl);
     return($matchV, \$sStri, undef) if($matchV >0);
   }
   return($mismV, undef, undef);
 }
 
 
 
 
 sub getMaxPairwiseMatch {
   
   my($curPloci, $curLoci, $penaltyFl, $penaltyFa
	, $queryStartI, $targetStartI ) = @_;
   
   $queryStartI=6 if(! defined($queryStartI));
   $targetStartI=8 if(! defined($targetStartI));

   
   
   my $cTotal=0;
   my $cSuccess=0;
   my $maxCEV=0;
   my $cOri=1;
   my $bestHit=undef;
   
   my @lkey= keys %{$curLoci};
   (@lkey == 1) or die "Error: There are more than one Hash-Entries for a queryEntry (getMaxPairwiseMatch)\n";
   my @plkey= keys %{$curPloci};
   (@plkey == 1) or die "Error: There are more than one Hash-Entries for a targetEntry (getMaxPairwiseMatch)\n";
   
   foreach my $ctE (@{$curPloci->{$plkey[0]}}) {
     # go through all members of the current Entry and keep only the best pairwise match with query item
     my $ctE_sDNA = $ctE->[$targetStartI];
     my $ctE_eDNA = $ctE->[$targetStartI +1];
     my $ctE_sAA = $ctE->[$queryStartI];
     my $ctE_eAA = $ctE->[$queryStartI +1];
     my $ctE_ori = -1;
     $ctE_ori = 1 if($ctE_sAA < $ctE_eAA);
     
		#match?
		foreach my $cqE (@{$curLoci->{$lkey[0]}}) {
		    my $cqE_sDNA = $cqE->[$targetStartI];
		    my $cqE_eDNA = $cqE->[$targetStartI +1];
		    my $cqE_sAA = $cqE->[$queryStartI];
		    my $cqE_eAA = $cqE->[$queryStartI +1];
		    my $cqE_ori = -1;
		    $cqE_ori = 1 if($cqE_sDNA < $cqE_eDNA);
		    
		    ($cqE_sAA, $cqE_eAA, $cqE_sDNA, $cqE_eDNA) = ($cqE_eAA, $cqE_sAA, $cqE_eDNA, $cqE_sDNA) unless($cqE_ori > 0);
		    
		    my $cE_matchV=0;
		    
		    if( ($ctE_sDNA == $cqE_sDNA) && ($ctE_eDNA == $cqE_eDNA) && ($ctE_sAA == $cqE_sAA) && ($ctE_eAA == $cqE_eAA)) {
			### pairwise match in positions???
			++$cTotal;
			
			if(($penaltyFl) && ($ctE_ori != $cqE_ori)) {
			    $cE_matchV += (($cqE->[-1]) * $penaltyFa / 100.0);
			} else {
			    $cE_matchV += ($cqE->[-1]); 
			    ++${cSuccess};
			}
			
			# check if current matchV is best
			if(${cE_matchV} > $maxCEV) {
			    $maxCEV= ${cE_matchV};
			    $bestHit = $ctE;
			    $cOri= $ctE_ori;
			    #$mCEVtotal=${cE_total};
			    #$mCEVsucc=${cE_succ};
			}
		    }
		}
	    }
	    
	    return($maxCEV, $bestHit, $cSuccess, $cTotal, $cOri);
	}
      

sub filterOrientation {
    my($liref, $orie
	,  $queryStartI, $targetStartI) = @_;
	
   $queryStartI=6 if(! defined($queryStartI));
   $targetStartI=8 if(! defined($targetStartI));

    my $resref;
    @{$resref}=();
    for(my $lii=0; $lii<@{$liref}; ++$lii) {
	my $curEntry = $liref->[$lii];
	my $curOri = 0;
	$curOri = 1 if(SynBlast::MyUtils::getSgn($curEntry->[$targetStartI +1] - $curEntry->[$targetStartI]) == SynBlast::MyUtils::getSgn($curEntry->[$queryStartI +1] - $curEntry->[$queryStartI]));
	if($orie == $curOri) {
	    push(@{$resref}, $curEntry);
	}
    }
    return($resref);
}

	
 
 sub buildBlastHitSet {
    
     # liref is supposed to be sorted descending and according to end-positions, with idx sortIdx
     my($liref, $sortIdx, $mOvl, $maxSize, $meanFl
		,  $verboseFl
		,  $queryStartI, $targetStartI) = @_;

    $queryStartI=6 if(! defined($queryStartI));
    $targetStartI=8 if(! defined($targetStartI));
    
    $meanFl = 0 if(! defined($meanFl));

         my @targetRegions=();
	 my $strtIdx=0;
	 while($strtIdx < @{$liref}) {
	     ### start new set of regions with current Ori (as distance on dna-level was to large)
	     my $curEntry = $liref->[$strtIdx];
	  	 
	     my $lastRend=$curEntry->[$sortIdx];
	     my $lastRstart=$curEntry->[$sortIdx -1];
	     my @curReg=();
	     my $groupKey = ( $curEntry->[$sortIdx-1] . "_" . $curEntry->[$sortIdx] );
	     if($meanFl) {
		 my $meanPos=$curEntry->[$sortIdx-1] - 1 + ($curEntry->[$sortIdx] - $curEntry->[$sortIdx-1] + 1)/2.;
		 $groupKey = ( $meanPos . "_" . $meanPos );
	     }
	     push(@curReg, { ( $groupKey => [ ($curEntry) ]) } );
	     push(@targetRegions, \@curReg);
	     
	     
	     my $li=$strtIdx+1; 
	     while($li< @{$liref}) {
		 my $cEi= $liref->[$li];
	
		my $frac_to_ext=0.5; #0.3  # specifies what fraction of maxSize is allowed at maximum as distance between neighbored hits triggering an extension of region size
		### region size in some cases has to be larger in order to get the relevant hits into one group
		 my $potentialNewSize=($lastRend - $cEi->[$sortIdx -1]);
		 my $curDista=$lastRstart - $cEi->[$sortIdx] ; 
		 if( ($sortIdx == $targetStartI+1)
			&& ($potentialNewSize >= $maxSize)  ) {###  potential new region size exceeds maxSize?
				
				#
				##print STDERR "exceeding maxSize ${maxSize}: ${potentialNewSize}\n";
				#
				
			if($curDista <= int($maxSize * $frac_to_ext))   {### but distance to neighbored hit is at max $frac_to_ext*maxSize -> then do extend to larger size
				warn "blasthit region size extended to ${potentialNewSize}, because dist (${curDista}) <= ${frac_to_ext}*maxSize (" 
					.  int($maxSize* $frac_to_ext) . ")!!!\n"
					. join("  ", @{$cEi}) . "\n" if($verboseFl);
			} else {  #distance to neighborhood too "large"?
				$strtIdx=$li;	### make up a new targetRegion in outer loop
				last;
			}
		}
		     
		 my $groupDone=0;
		 ### group entry in current last targetRegion in List
		 for(my $cli=0; $cli< @{$targetRegions[-1]}; ++$cli) {
		     
		     ## only last entry to look at - maybe sufficient
		     
		     #foreach my $key (keys %{$targetRegions[-1]}) {
		     my @kes = keys %{$targetRegions[-1]->[$cli]};
		     die "Error: only one key allowed per array here!\n" unless(@kes == 1);
		     my $key=$kes[0];
		     my @key = split(/_/, $key);  ### dnaStart_dna_End
		     
		     $groupKey = ( $cEi->[$sortIdx-1] . "_" . $cEi->[$sortIdx] );
		     my $meanPos = undef;
		     if($meanFl) {
			 $meanPos=$cEi->[$sortIdx-1] - 1 + ($cEi->[$sortIdx] - $cEi->[$sortIdx-1] + 1)/2.;
			 $groupKey = ( $meanPos . "_" . $meanPos );
		     }
		     
		     ## check next entry if overlap not larger than allowed, i.e. if there is no need to combine/group this entry to current one
		     my $curOvlS= SynBlast::MyUtils::getOverlapSize($key[0], $key[1], $cEi->[$sortIdx -1], $cEi->[$sortIdx]);
		     next if((! $meanFl) 
				&& ( $curOvlS <= $mOvl) 
				&& ( abs($cEi->[$sortIdx -1] - $cEi->[$sortIdx]) + 1 >=  2 * $curOvlS )  ## new interval is in limit of overlap 
														#but also longer in size as actual_overlap (at least twice as long)
				);
		
		     next if( $meanFl && ($key[0] > $meanPos)); ## case if meanPos larger leads to grouping even though the meanpos is different?
		     ## else
		     ## calculate new boundaries (only dna_end should be changed!)
		     if(! $meanFl) {
			 my $newKey = SynBlast::MyUtils::min($key[0], $cEi->[$sortIdx -1]) . "_" . SynBlast::MyUtils::max($key[1], $cEi->[$sortIdx]);
			 print STDERR "WARNING/ERROR: old regionEnd was " . $key[1] . ", new is " . $cEi->[$sortIdx] . "!\n" if($key[1] < $cEi->[$sortIdx]);
			 ##insert
			 if($newKey ne $key) {
			     my ($oldValue) = delete($targetRegions[-1]->[$cli]->{$key});
			     $targetRegions[-1]->[$cli]->{$newKey} = $oldValue;
			     #   print STDERR "Change in regKey: ${key} vs ${newKey}\n";
			     #  print STDERR ";";
			 } else {
			     #  print STDERR "no change in regKey...\n";
			     #print STDERR ".";
			 }
			 push(@{$targetRegions[-1]->[$cli]->{$newKey}}, $cEi);
		 
		     } else { ## meanpos is the same or larger?
			 print "WARNING: meanpos was larger (${meanPos} vs. " . $key[0] . ")\n" if($meanPos > $key[0]);
			 my $newKey = $key;
			 push(@{$targetRegions[-1]->[$cli]->{$newKey}}, $cEi);
		     }
		     ++$groupDone;
		     last; ## stop here as current entry inserted...
		 }
		 
		 
		 if(! $groupDone) {
		     ### not yet inserted?
		     $groupKey = ( $cEi->[$sortIdx-1] . "_" . $cEi->[$sortIdx] );
		     my $meanPos = undef;
		     if($meanFl) {
			 $meanPos=$cEi->[$sortIdx-1] - 1 + ($cEi->[$sortIdx] - $cEi->[$sortIdx-1] + 1)/2.;
			 $groupKey = ( $meanPos . "_" . $meanPos );
		     }
		     my $newKey = $groupKey;
		     push(@{$targetRegions[-1]}, { ($newKey => [ ($cEi) ] ) } );
		     $lastRstart=SynBlast::MyUtils::min($lastRstart, $cEi->[$sortIdx -1]);
		     ++$groupDone;
		 }
		 ###
		 ++$li;
	     }
	     $strtIdx = $li if($li >= @{$liref});
	 }
     
     return(\@targetRegions);
 }



}



##########################################################################
# HSP chaining procedure old version (deep recursions are here possible!)
##########################################################################
sub calcSplAln {
    my ($lref, $resultRef, $queryID, $maxSize, $scutoff, $maxOvl, $verboseFl, $queryL
       , $queryStartI, $targetStartI, $lenI, $preLcds, $posDiffFac) = @_;
       
    # maxSize is maximal length of range for loci (should depend upon reference-size)
    # scutoff is the cutoff for the bitscore necessary for any valid putative final loci
    # maxOvl is the max. overlap length (protein overlap) allowed in putative locis (consisting of several (overlapping) entries/blast-hits
    $verboseFl=0 if(! defined($verboseFl));

$posDiffFac=3 unless(defined($posDiffFac));

  $queryStartI=6 if(! defined($queryStartI));
  $targetStartI=8 if(! defined($targetStartI));
  $lenI=3 if(! defined($lenI));
  

# assumption!!!: lref already sorted ascending acc. to target/DNA-end positions !!!
    my $dummy=0;
    for(my $li=0; $li< @{$lref}; ++$li) {
	if($lref->[$li]->[$targetStartI+1] >= $dummy) {
	    $dummy=$lref->[$li]->[$targetStartI+1];
	} else {
	    die "ERROR: lref not sorted asc. in target-end-position!\n";
	}
    }




    my $scMax=$scutoff;
    
    my @todoList=@{$lref};
    while((@todoList > 0) && ($scMax >= $scutoff)) {
	my @cList= @todoList;
	my $todoN=@todoList; # number of entries to do
	$scMax=0;
	my @selEntries=();
	my $finalOri=0;
	
	for(my $jj=$todoN -1; $jj >= 0; --$jj) {
	    my $lastRef=pop(@cList);
	    my @selE=();
	
	    my ($cScore, $cOri) =myScFun($lastRef, \@cList, \@selE, $lastRef->[$targetStartI+1], $maxSize, $maxOvl , $queryStartI, $targetStartI, $posDiffFac);#, \@nselE);
	    next if($cScore < $scutoff);
	    if($scMax <= $cScore) {
		$scMax = $cScore;
		@selEntries=@selE;
		$finalOri=$cOri;
	    }
	}
	last if($scMax < $scutoff);
	cleanUptdl(\@todoList, \@selEntries); # removes all entries in final-selection out of todoList (~difference)	
	print "score($queryID)=" . $scMax . "\n" if($verboseFl);
	SynBlast::MyUtils::printListList(\@selEntries, "\t", "\n") if($verboseFl > 1);
	print "---\n" if($verboseFl > 1);
	SynBlast::MyUtils::printListList(\@todoList, "\t", "\n") if($verboseFl > 1);
	print "ERROR: there should be " . ($todoN - @selEntries) . " in todo-list!\n\n" if(@todoList != (($todoN - @selEntries)));

	my @curRangeField=getRangeField(\@selEntries, $queryStartI, $targetStartI, $lenI, $preLcds, $posDiffFac);
	my %hcontent=();
	$hcontent{'score'}=$scMax;
	$hcontent{'hits'}=\@selEntries;
	$hcontent{'ori'}=$finalOri;
	$hcontent{'startPos'}=$curRangeField[0];
	$hcontent{'endPos'}=$curRangeField[1];
	$hcontent{'hitsCount'}= $curRangeField[2];
	$hcontent{'meanPos'}=$curRangeField[3];
	$hcontent{'hitsLength'}= $curRangeField[4];
	$hcontent{'queryCov'} = ( int($queryL) ? ("" . (int(( $curRangeField[5] / $queryL * 1000) + 0.5) / 10) . "%") : ("n.a."));
      

	$hcontent{'chr'}= $curRangeField[6];
	$hcontent{'meanIdentity'}= ( "" . $curRangeField[7] ."%"); 
	$hcontent{'CDSannot'} = $curRangeField[8];
	$hcontent{'queryID'} = $selEntries[0]->[0];
	$hcontent{'uncoveredAA1'} = ( int($queryL) ?  ( ($curRangeField[9]->[0] -1) . "(" . (int( ( ($curRangeField[9]->[0] -1) / $queryL * 1000.) + 0.5) /10.) . "%)" ) : ( ($curRangeField[9]->[0] -1) )); #  ("n.a."));
	$hcontent{'uncoveredAA2'} = ( int($queryL) ?  ( ($queryL - $curRangeField[9]->[1]) . "(" . (int( ( ($queryL - $curRangeField[9]->[1]) / $queryL * 1000.) + 0.5) /10.) . "%)" ) : ("n.a.")); #( ($queryL - $curRangeField[9]->[1]) . "aa")); # ("n.a."));
	## put here check on overlap with prior locis for queryID
	if(checkLociOverlapOK(\@{$resultRef->{$queryID}}, \@curRangeField, $scMax)) {
	    push(@{$resultRef->{$queryID}}, \%hcontent);
	} else {
	    print "Putative Loci \"" . join("_", @curRangeField) . "\t${scMax}\" is overlapping with a higher scoring loci -> skipping this entry...\n" if($verboseFl);
	}
    }
    print "END\n\n" if($verboseFl);
    

    sub cleanUptdl {
     my ($todoref, $selref) = @_;
     my @IdxToRemove=();
     for(my $jj=0; $jj<@$todoref; $jj++) { # IdxToRemove is filled in ascending order!
	 foreach my $curE (@$selref) {
	     push(@IdxToRemove, $jj) if(join("_", @$curE) eq join("_", @{$todoref->[$jj]}));
             # curE is found in todoref at pos jj
	 }
     }
     (@IdxToRemove <= @$todoref) or die "Error: more idxs to remove than available...\n";
     for(my $ii=@IdxToRemove -1; $ii>=0; $ii--) {
	 #note: descending order necessary as upper idxpositions in todoref changing after splice...
	 splice(@$todoref, $IdxToRemove[$ii], 1); 
     }
 }

    
    sub myScFun {
	my ($cHit, $lref, $selref, $lociStartEnd, $maxSize, $maxOvl
	, $queryStartI, $targetStartI, $posDiffFac) = @_;
	# lociStartEnd is to calculate the rangeSize (largest target-endpos in current range)
	
	$posDiffFac=3 unless(defined($posDiffFac));
	
	$queryStartI=6 if(! defined($queryStartI));
        $targetStartI=8 if(! defined($targetStartI));

	my $cTstart=$cHit->[$targetStartI];
	my $cTend=$cHit->[$targetStartI+1];
	my $cQstart=$cHit->[$queryStartI];
	my $cQend=$cHit->[$queryStartI+1];
	my $cQori=-1;
	$cQori=1 if ($cQstart <= $cQend);
	my $cMaxFunRes=0;
	for(my $ii=@{$lref} -1; $ii>=0; --$ii) {
	    my $cH= $lref->[$ii];
	    my $cOverlapP=0;
	    $cOverlapP= $cQstart - $cH->[$queryStartI+1] + 1 if(($cH->[$queryStartI] > $cH->[$queryStartI+1]) && ($cQori < 0));
	    $cOverlapP= $cH->[$queryStartI+1] - $cQstart + 1 if(($cH->[$queryStartI] <= $cH->[$queryStartI+1]) && ($cQori > 0));
	    $cOverlapP=0 if ($cOverlapP < 1);
	   # if( ($cH->[9] < $cTstart) && 
	    if( ($cOverlapP <= $maxOvl) && ($cH->[$targetStartI+1] < $cTstart + ($posDiffFac * $cOverlapP)) &&
		( ( ($cH->[$queryStartI+1] - $cOverlapP < $cQstart) && ($cH->[$queryStartI] <= $cH->[$queryStartI+1]) && ($cQori > 0) ) 
		  || ( ($cH->[$queryStartI+1] + $cOverlapP > $cQstart) && ($cH->[$queryStartI] > $cH->[$queryStartI+1]) && ($cQori < 0) )
		  ) 
		) {
		#print "# cH is of desired property (chain-like)...\n";
		
		#check if size of curLoci is in limit...
		if(abs($lociStartEnd - $cH->[$targetStartI]) + 1 <= $maxSize) {
		    
		    
		    #calculate myScFun
		    my @clist= @$lref;
		    my @subsel=();
		    #my @nsubsel=();
		    splice(@clist, $ii);#, 1); 
		    my ($cSc, $subOri) = myScFun($cH, \@clist, \@subsel, $lociStartEnd, $maxSize, $maxOvl , $queryStartI, $targetStartI, $posDiffFac);
		    if($cMaxFunRes <= $cSc) {
			$cMaxFunRes = $cSc;
			@{$selref} = @subsel;
			die "Error (myScFun): match-Orientation in chain is not consistent!\n" unless ($subOri eq $cQori);
		    }
		}
	    }
	}
	# maximum subselection down from cHit calculated... add cHit to selection and return final score...
	push(@$selref, $cHit);
	return((($cHit->[-1]) + $cMaxFunRes, $cQori));
    }
}




sub checkLociOverlapOK {
       my ($listref, $rangeField, $score) = @_;
       my $ovlCount=0;
	for(my $li=0; $li<@{$listref}; $li++) {
	    my $curHaRef= $listref->[$li];
	    if( SynBlast::MyUtils::getOverlapSize($rangeField->[0], $rangeField->[1], $curHaRef->{'startPos'}, $curHaRef->{'endPos'}) != 0) {
		if(defined($score) && ($score > $curHaRef->{'score'})) {
		    warn "current Loci is overlapping and has a better score!!!!!!! Overlapping existing entry will be deleted...\n";
		    $listref->[$li] = undef;
		} else {
		    ++$ovlCount;
		}
	    }
	}
	return(($ovlCount == 0));
    }
    



sub getRangeField {
    my ($selref
	, $queryStartI, $targetStartI, $lenI, $preLcds, $posDiffFac) = @_;
	
	$posDiffFac=3 unless(defined($posDiffFac));

	#$preLcds= 3*300 unless(defined($preLcds));
	$preLcds= $posDiffFac * 300 unless(defined($preLcds));
	
	$queryStartI=6 if(! defined($queryStartI));
	$targetStartI=8 if(! defined($targetStartI));
	
	$lenI=3 if(! defined($lenI)); ## aln.Length index position
	
	(@$selref > 0) or return( ("nix", "nix", 0, -1, 0, 0, "noChromosome", 0) );
	
	my $rstart=$selref->[0]->[$targetStartI];
	my $rend=0;
	my $meanP=0;
	my %sortH=();
	my $chrName= $selref->[0]->[1];
	
	
	for(my $i=0; $i<@{$selref}; $i++) {
	    $rstart=SynBlast::MyUtils::min($rstart, $selref->[$i]->[$targetStartI]);
	    $rend=SynBlast::MyUtils::max($rend, $selref->[$i]->[$targetStartI+1]);
	    push( @{$sortH{$selref->[$i]->[$targetStartI]} } , $selref->[$i]);
	}
	
	my @sortBlHits=();
	my $totalMatchl=0;

    my $totProtL = 0; # for the number of AAs covered by BlastHit
    my @rangeAA=(undef, undef); # here the range of protein covered; to calculate later the number of AAs out of covered range up and downstream of covered protein range
    
    my @lastProtPos=(0,0);
    my $totIdL = 0; ## counts the number of identical matches
    my $totAlnL = 0; ## sums up the aln length column (3) of the entries
    my @CDSannot=(); ## stores information for CDS positions (not necessarily true/single exons) (only check for overlap with previous CDSannotation's queryEnd to correct for overlaps
			## format is: [0]-> (targetStart, targetEnd), [1]->(queryStart, queryEnd), 
    
    
	my $negori=undef;
	foreach my $startPos (sort {$a <=> $b} (keys %sortH)) {
	    foreach my $liref (@{ $sortH{$startPos} }) {
		push(@sortBlHits, [ @{$liref} ]);
		$totalMatchl+= (abs($liref->[$targetStartI+1] - $liref->[$targetStartI]) + 1);
		
		#my $negori=0;
		if(! defined($negori)) {
			$negori=0;
			$negori=1 if($liref->[$queryStartI] > $liref->[$queryStartI+1]);
		}
		die "ERROR: different orientation within HSP chain!\n" if((! $negori) && ($liref->[$queryStartI] > $liref->[$queryStartI+1]));
		die "ERROR: different orientation within HSP chain!\n" if(($negori) && ($liref->[$queryStartI] < $liref->[$queryStartI+1]));

		
		
		my $zuwachs= abs($liref->[$queryStartI] - $liref->[$queryStartI+1]) + 1 - SynBlast::MyUtils::getOverlapSize(@lastProtPos,  $liref->[$queryStartI], $liref->[$queryStartI+1]);
		$zuwachs= (abs($liref->[$queryStartI] - $liref->[$queryStartI+1]) + 1 - SynBlast::MyUtils::getOverlapSize(@lastProtPos,  $liref->[$queryStartI+1], $liref->[$queryStartI])) if($negori);  ## orientation -;
		$totProtL += $zuwachs;
		@lastProtPos=($liref->[$queryStartI], $liref->[$queryStartI+1]);
		@lastProtPos=($liref->[$queryStartI+1], $liref->[$queryStartI]) if($negori);  ## orientation -
		
		@rangeAA= ( SynBlast::MyUtils::min($rangeAA[0], $lastProtPos[0]), SynBlast::MyUtils::max($rangeAA[1], $lastProtPos[1]) );
		
		### idx2= perc.identity, idx3=aln.length
		$totIdL += int(($liref->[2] * $liref->[$lenI] / 100 ) + 0.5);
		$totAlnL += int($liref->[$lenI]);
		my $overlapAA=0;
		if(@CDSannot) {
				$overlapAA= $CDSannot[-1]->[2]->[1] - $liref->[$queryStartI] +1;
				# last CDS queryEndPosition - currentCDS queryStartPos +1
				$overlapAA= $liref->[$queryStartI] - $CDSannot[-1]->[2]->[1] + 1 if($negori);
			}
		my $missCh="";
		if($overlapAA < 0) { ## no real overlap, but space inbetween (no contiguous mapping)
			$overlapAA=0;
			$missCh="<";
		}
		if($negori) {
			push(@CDSannot, [ ("${missCh}" ,  [ ( ($liref->[$targetStartI] + $posDiffFac*$overlapAA), $liref->[$targetStartI+1]) ], [ ($liref->[$queryStartI] - $overlapAA, $liref->[$queryStartI+1]) ] ) ] );
		} else {
			push(@CDSannot, [ ( "${missCh}" , [ ( ($liref->[$targetStartI] + $posDiffFac*$overlapAA), $liref->[$targetStartI+1]) ], [ ($liref->[$queryStartI] + $overlapAA, $liref->[$queryStartI+1]) ] ) ] );
		}
				
				
	    }
	}
	
		
	#my $totalNt= $CDSannot[-1]->[1]->[1]-$CDSannot[0]->[1]->[0] + 1 + 2*$preLcds;
	my @totalInterv= ($CDSannot[0]->[1]->[0] - $preLcds, $CDSannot[-1]->[1]->[1] + $preLcds);
	my $expectLen=$totalInterv[1]-$totalInterv[0]+1;
	$totalInterv[0] = 1 if($totalInterv[0] < 1); # if startPos is smaller than preLcds, use only available pre-sequence (from start of contig)
	push(@totalInterv, $totalInterv[1]-$totalInterv[0]+1); # total intervall length
	push(@totalInterv, $expectLen); # expected intervall length
	
	
	my %relCDSannot=();
	my $origCDSoriNeg=0;
	$origCDSoriNeg=1 if($CDSannot[0]->[2]->[0] > $CDSannot[0]->[2]->[1]);
	if($origCDSoriNeg) {
		for(my $lli=@CDSannot-1; $lli>=0; --$lli) {
			push(@{$relCDSannot{"rel"}}, [    ( $CDSannot[$lli]->[0] ? ">" : "")
					, [   ($totalInterv[2] - ($CDSannot[$lli]->[1]->[1]- $totalInterv[0]+1) +1) ### startPos rev.complement
					    , ($totalInterv[2] - ($CDSannot[$lli]->[1]->[0]- $totalInterv[0]+1) +1) ## endPos rev.complement
					  ]
					, [ reverse @{$CDSannot[$lli]->[2]} ]
				] );
		}
		$relCDSannot{"revFl"} = 1;
	} else {
		foreach my $cAnnot (@CDSannot) {
			push(@{$relCDSannot{"rel"}}, [    ( 
										$cAnnot->[0]
										, [ (  $cAnnot->[1]->[0] - $totalInterv[0]+1
											, $cAnnot->[1]->[1] - $totalInterv[0]+1
											)
											]
										, [ ( @{$cAnnot->[2]} )  ]
					)
				   ] );
		}
		$relCDSannot{"revFl"} = 0;
	}
	
	$relCDSannot{"orig"} =  \@CDSannot;  ### besser kopieren!!!
	$relCDSannot{"srcInterv"} = \@totalInterv;
	

	my $foundFl=0;
        my $preDNAl = int($totalMatchl /2) - 1;

    
	for(my $li=0; $li<@sortBlHits; $li++) {

            if(! $foundFl) {
		my $matchLength= abs($sortBlHits[$li]->[$targetStartI+1] - $sortBlHits[$li]->[$targetStartI]) + 1;
		if( $preDNAl > $matchLength) {
		    $preDNAl-=$matchLength;
		} else {
		    $foundFl++;
		    $meanP= $sortBlHits[$li]->[$targetStartI] + $preDNAl;
		}
	    }
	}
	
	return( ($rstart
		, $rend
		, int(@{$selref})
		, $meanP
		, $totalMatchl
		, $totProtL
		, $chrName
		, (int(($totIdL / $totAlnL * 1000) + 0.5) / 10.0)
		, \%relCDSannot
		, \@rangeAA
		)
	    );
	  #### mean percent identity containes error, as alns are partly overlapping and position of matches not known, hence matches in overlapping parts are count twice sometimes
	  #### therefore dividing the total Number of Matches through the total alignment length (including overlapping) 
    }
    







__END__

=head1 SynBlast::MySyntenyAln

SynBlast::MySyntenyAln::getSyntenyAln
SynBlast::MySyntenyAln::getBestBlastLoci
SynBlast::MySyntenyAln::getCombPlociRef
SynBlast::MySyntenyAln::calcSplAln2
SynBlast::MySyntenyAln::calcSplAln


=head1 SYNOPSIS

 require GD;
 
=head1 DESCRIPTION

functions implementing the synteny region evaluation/alignment part of the pipeline

=head1 AUTHOR

Joerg Lehmann, <lt>Ejoe@bioinf.uni-leipzig.deE<gt>

=cut

