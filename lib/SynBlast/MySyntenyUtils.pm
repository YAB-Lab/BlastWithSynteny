#!/usr/bin/env perl
# Last changed Time-stamp: <2008-08-19 15:06:06 joe>
# joe@bioinf.uni-leipzig.de
########################################################################

use strict;
use warnings;
use SynBlast::MyUtils;

package SynBlast::MySyntenyUtils; 
require Exporter;
my @ISA = qw(Exporter);
my @EXPORT = qw(blastfile2Lists infofile2Exonlist info2Infofile blastfile2sortedLists getSynRegInfo);
my @EXPORT_OK = @EXPORT;



sub getSynRegInfo {
my ($blastFile) = @_;
open(BLASTF, '<', "${blastFile}")  or die "Error while reading ${blastFile}: $!\n";
my $counter=0;
my %infoH=();
while(<BLASTF>) {
    if ($_ =~ m/^\#/) {
	chomp;
	my @curF=split(/\s+/);
	if(@curF > 1) {
	    substr($curF[0], 0, 1) = "" if(substr($curF[0], 0, 1) eq "#");
	    $infoH{$curF[-2]}=$curF[-1];
	}
    }
}

return($infoH{'ensemblDate'},
       $infoH{'synOrg'},
       $infoH{'synChr'},
       $infoH{'synStart'},
       $infoH{'synEnd'});
}      







# reads a blastfile as-is into hash containing a list of entries grouped by idx $keyIdx
sub blastfile2Lists {
    my ($blastFile, $hashref, $keyIdx) = @_;
    defined($keyIdx) or die "Error: no grouping index selected (blastfile2Lists)\n";
    open(BLASTF, '<', "${blastFile}")  or die "Error while reading ${blastFile}: $!\n";
    my $counter=0;
    while(<BLASTF>) {
	if ($_ =~ m/^[^\#]/) {
	    chomp;
	    my @curF=split(/\s+/);
	    ($keyIdx < @curF) or die "Error: grouping index out of range (blastfile2Lists)\n";
	    push(@{$hashref->{$curF[$keyIdx]}}, [ @curF ]);
	    $counter++;
	}
    }
    close(BLASTF);
    return($counter);
}





# reads a blastfile as-is into hash containing a SORTED list of entries (hash grouped by idx $keyIdx; 
# each list sorted according to $sortIdx
# optional filtering: allowing only entries with $fvalue at $findex
sub blastfile2sortedLists {
    my ($blastFile, $hashref, $keyIdx, $sortIdx, $findex, $fvalue, $swapFl, $keyIdx2, $queryIdx, $targetIdx, $wuformatFl) = @_;
    
    defined($keyIdx) or die "Error: no grouping index selected (blastfile2sortedLists)\n";
    defined($sortIdx) or die "Error: no sort index selected (blastfile2sortedLists)\n";
    (! defined($findex)) or (defined($fvalue)) or die "Error: filter-value not defined! (blastfile2sortedLists)\n";
    if(! defined($swapFl)) { $swapFl=1; }
    
    my @blastFiles;
    if(ref $blastFile) {
      @blastFiles = @{$blastFile};
    } else {
      @blastFiles = ($blastFile);
    }
    
    $queryIdx = 6 if(! defined($queryIdx));
    $targetIdx = 8 if(! defined($targetIdx));
    
    if($wuformatFl) {
	    $queryIdx = 8;
	    $targetIdx = 10;
    }
    
    
    my $counter=0;
    my %blastsortH=();
    my $dum=0;
    #$dum=2 if(! $swapFl);
    $dum=$targetIdx-$queryIdx if(! $swapFl);
    
    foreach $blastFile (@blastFiles) {
      open(BLASTF, '<', "${blastFile}")  or die "Error while reading ${blastFile}: $!\n";
      
      while(<BLASTF>) {
	if ($_ =~ m/^[^\#]/) {
	  chomp;
	  my @curF=split; #(/\s+/);
	  next if(! @curF);
	  if($wuformatFl) {  ### 22 strd fields, optionally two more (24)
		  next unless((@curF >= 22) && (@curF <= 24));
		  # Fields: qid   sid     E       N       Sprime  S       alignlen        nident  npos    nmism   pcident pcpos   qgaps   qgaplen sgaps   sgaplen qframe  qstart  qend    sframe  sstart  send    group   links
		   @curF = ( $curF[0]  # pos0: queryID
					,$curF[1] #pos1: subjectID
					,$curF[10]  #pos2:perc.id
					,$curF[11] #pos3: perc.pos
					,$curF[16] . "/" . $curF[19] #pos4: qry/sbj-frame
					,$curF[6]  #pos5: aln-length
					,$curF[9]   #pos6: nmism
					,$curF[14]  #pos7: gap opens?
					,$curF[17]  #pos8:q.start
					,$curF[18]   #pos9:q.end
					,$curF[20]   #pos10:s.start
					,$curF[21]   #pos11:s.end
					,$curF[2]    #pos12:evalue
					,$curF[4]     #pos13:bitscore
				);
	  }
	  
	  
	  
	 die "ERROR: Could not read blast file properly! NCBI tabular blast output format expected (14 columns)!\n" unless(@curF < 15 );
	  		  
	  ($sortIdx < @curF) or die "Error: sort index out of range (blastfile2sortedLists; ${sortIdx} should be less than " . @curF . ")\n";
	  (! defined($findex)) or ($findex < @curF) or die "Error: filter-index out of range (blastfile2sortedLists)\n";	    
	  
	  do { SynBlast::MyUtils::swapArray(\@curF, $queryIdx, $queryIdx+1);
	    SynBlast::MyUtils::swapArray(\@curF, $targetIdx, $targetIdx+1);
	    #my $tmp=$curHit[6]; $curHit[6]=$curHit[7]; $curHit[7]=$tmp; $tmp=$curHit[8]; $curHit[8]=$curHit[9]; $curHit[9]=$tmp;
	  } unless ( int($curF[$targetIdx+1 -$dum]) >= int($curF[$targetIdx -$dum]));
	  
	  
	  if((! defined($findex)) || ($curF[$findex] eq "${fvalue}")) {
	    push(@{$blastsortH{$curF[$sortIdx]}}, [ @curF ]);
	  }
	}
      }
      close(BLASTF);
    }
    
    ### sort according to sortIdx
    my @sortKeys= sort {$a <=> $b} (keys %blastsortH);
    for(my $li=0; $li < @sortKeys; ++$li) {
      foreach my $cEref (@{$blastsortH{$sortKeys[$li]}}) {
	($keyIdx < @{$cEref}) or die "Error: grouping index out of range (blastfile2sortedLists)\n";
	if(defined($keyIdx2)) {
	  # a second keyIdx to group on?
	  ($keyIdx2 < @{$cEref}) or die "Error: grouping index2 out of range (blastfile2sortedLists)\n";
	  push(@{$hashref->{$cEref->[$keyIdx]}->{$cEref->[$keyIdx2]}}, $cEref);
	} else {
	  push(@{$hashref->{$cEref->[$keyIdx]}}, $cEref);
	}
	$counter++;
      }
    }
    
    return($counter);
  }






sub lists2sortedLists {
    my ($inpHash, $hashref, $keyIdx, $sortIdx, $findex, $fvalue, $swapFl
        , $keyIdx2, $queryIdx, $targetIdx) = @_;
	
    defined($keyIdx) or die "Error: no grouping index selected (lists2sortedLists)\n";
    defined($sortIdx) or die "Error: no sort index selected (lists2sortedLists)\n";
    (! defined($findex)) or (defined($fvalue)) or die "Error: filter-value not defined! (lists2sortedLists)\n";
    if(! defined($swapFl)) { $swapFl=1; }
    
    
    $queryIdx = 6 if(! defined($queryIdx));
    $targetIdx = 8 if(! defined($targetIdx));
    
    my $counter=0;
    my %blastsortH=();
    my $dum=0;
    #$dum=2 if(! $swapFl);
    $dum=$targetIdx-$queryIdx if(! $swapFl);
    
    foreach my $inpKey (keys %{$inpHash}) {
	foreach my $curF (@{$inpHash->{$inpKey}}) {
	    my @curF=@{$curF};
	    ($sortIdx < @curF) or die "Error: sort index out of range (lists2sortedLists)\n";
	    (! defined($findex)) or ($findex < @curF) or die "Error: filter-index out of range (lists2sortedLists)\n";	    
	    
	    do { SynBlast::MyUtils::swapArray(\@curF, $queryIdx, $queryIdx+1);
	      SynBlast::MyUtils::swapArray(\@curF, $targetIdx, $targetIdx+1);
	  } unless ( int($curF[$targetIdx+1-$dum]) >= int($curF[$targetIdx-$dum]));

	    
	    if((! defined($findex)) || ($curF[$findex] eq "${fvalue}")) {
		push(@{$blastsortH{$curF[$sortIdx]}}, [ @curF ]);
	    }
	}
    }
    ### sort according to sortIdx
    my @sortKeys= sort {$a <=> $b} (keys %blastsortH);
    for(my $li=0; $li < @sortKeys; ++$li) {
	foreach my $cEref (@{$blastsortH{$sortKeys[$li]}}) {
	    ($keyIdx < @{$cEref}) or die "Error: grouping index out of range (lists2sortedLists)\n";
	    
	    if(defined($keyIdx2)) {
		# a second keyIdx to group on?
		($keyIdx2 < @{$cEref}) or die "Error: grouping index2 out of range (lists2sortedLists)\n";
		push(@{$hashref->{$cEref->[$keyIdx]}->{$cEref->[$keyIdx2]}}, $cEref);
		} else {
		push(@{$hashref->{$cEref->[$keyIdx]}}, $cEref);
	    }
	    $counter++;
	}
    }
    
    return($counter);
}












### creates a new InfoFile out of infoFile adding additional data columns out of addInfoHash (key is protID)
sub info2Infofile {
    my ($infoFileSrc
    , $extInfoF
    , $addInfoH
    , $verboseFl
    , $queryColName
    ) = @_;
    
    $queryColName = "protID" if(! defined($queryColName));
    
    
    ### infoFileSrc is used as source and is extended with info of addInfoHash (output to extInfoF)
#    * .protein.info  -> protein.infoW
    if(! defined($extInfoF)) {
#	$extInfoF = $infoFileSrc . "W";
	$extInfoF = $infoFileSrc
    }

    if(! (defined($addInfoH) && exists($addInfoH->{'dattype'})) ) {
	die "ERROR: No additional data available!\n";
    }
    
    $verboseFl=1 if(! defined($verboseFl));
    
    open(INFOFS, '<', "${infoFileSrc}")  or die "Error while reading ${infoFileSrc}: $!\n";
    my @newInfoCont=();

    my $headerl;
    my %infoH=();
    while(defined($headerl = <INFOFS>) && (${headerl} =~ m/^\#/)) {
    	if(${headerl} =~ m/^\#!/) {
    			chomp($headerl);
    			my @defl=split(/=/, $headerl);
    			$infoH{substr($defl[0], 2)} = $defl[-1];
    			push(@newInfoCont, $headerl . "\n") if($defl[0] ne "#!columns");
      } else {
      	#chomp($headerl);
      	push(@newInfoCont, $headerl)
      }
      
    }
    return(undef) if(! defined($headerl)); ## no content after description?
    foreach ( qw( columns focalgene organism seqregion from to )) {
    	die "Error: Could not find \"#!$_=\" definition in ${infoFileSrc}!\n" unless exists($infoH{$_});
    }
    
    #(defined($headerl = <INFOF>) && (${headerl} =~ m/^\#columns/)) or die "Error while reading ${infoFile}: 2nd line is expected to be the definition of the columns, starting with \"#columns\"\n";
    my @tcolumns = split(/,/, $infoH{"columns"});
    my %colIdx=();
    ### which idx/col is which info to find at?
    for(my $ti=0; $ti<@tcolumns; ++$ti) {
    		$colIdx{$tcolumns[$ti]}=$ti;
    }
    foreach ( $queryColName, "startPos", "endPos", "ori", "geneID", "protLength") {
    	die "Error: Could not find \"$_\" column in ${infoFileSrc}!\n" unless exists($colIdx{$_});
    }
    if(exists($colIdx{"selfScore"})) { ## selfScore already in infoFile...
    		push(@newInfoCont, "#!columns=" . join(",", @tcolumns) . "\n");
    	} else {
    		push(@newInfoCont, "#!columns=" . join(",", @tcolumns) . ",selfScore" . "\n");
    }
        
   my $colTotal=@tcolumns; #0;
   my $cPos=0;
   	do {
    	if ($headerl =~ m/^[^\#]/) {
  	    chomp($headerl);
	      my @cLi = split(/\s+/, $headerl);
	      die "ERROR: the number of data columns is not the same as in column-definition (pos${cPos})!\n" unless(@cLi eq $colTotal);
	      ++$cPos;
	 
	      if(exists($addInfoH->{$cLi[$colIdx{$queryColName}]})) {
					if(! exists($colIdx{"selfScore"})) { ## add new
						push(@newInfoCont, join("\t", @cLi, $addInfoH->{$cLi[$colIdx{$queryColName}]} ) . "\n");
		      } else { 
		    ### already data for selfComp available -> will be overwritten with new data,  if complete
		    warn("Note: already available selfComparison-score is updated/changed!\n") if($cLi[$colIdx{'selfScore'}] != $addInfoH->{$cLi[$colIdx{$queryColName}]});
		    $cLi[$colIdx{'selfScore'}] = $addInfoH->{$cLi[$colIdx{$queryColName}]};
		    push(@newInfoCont, join("\t", @cLi) . "\n");
		   }
	    } else {
		die "Could not find additional data (" . $addInfoH->{'dattype'} . ") for " . $cLi[$colIdx{$queryColName}] . "! InfoFile is not changed/created.\n";
    }
	} else {
	    push(@newInfoCont, $headerl);
	 }
 } while(defined($headerl=<INFOFS>));
    
    close(INFOFS);
    
    
    open(EINFOF, '>', "${extInfoF}")  or die "Error while creating ${extInfoF}: $!\n";
    print EINFOF join("", @newInfoCont); 
    close(EINFOF);
    print "${extInfoF} was successfully created/extended...\n" if($verboseFl);
    
    return(1);
}





sub infofile2Exonlist {
    my ($infoFile
      , $listref
      , $queryIdxHash
      , $refexonsList
      , $exonIDlist
      , $verboseFl
      , $queryColName) = @_;
    
    $verboseFl=1 if(! defined($verboseFl));
    $queryColName = "protID" if(! defined($queryColName));
    
    open(INFOF, '<', "${infoFile}")  or die "Error while reading ${infoFile}: $!\n";
  
    my $headerl;
    my %infoH=();
    while(defined($headerl = <INFOF>) && (${headerl} =~ m/^\#/)) {
    	if(${headerl} =~ m/^\#!/) {
    			chomp($headerl);
    			my @defl=split(/=/, $headerl);
    			$infoH{substr($defl[0], 2)} = $defl[-1];
      }
    }
    return(undef) if(! defined($headerl)); ## no content after description?
    foreach ( qw( columns focalgene organism seqregion from to )) {
    	die "Error: Could not find \"#!$_=\" definition in ${infoFile}!\n" unless exists($infoH{$_});
    }
    
    #(defined($headerl = <INFOF>) && (${headerl} =~ m/^\#columns/)) or die "Error while reading ${infoFile}: 2nd line is expected to be the definition of the columns, starting with \"#columns\"\n";
    my @tcolumns = split(/,/, $infoH{"columns"});
    my %colIdx=();
    ### which idx/col is which info to find at?
    for(my $ti=0; $ti<@tcolumns; ++$ti) {
    		$colIdx{$tcolumns[$ti]}=$ti;
    }
    foreach ( $queryColName, "startPos", "endPos", "ori", "geneID", "protLength") {
    	die "Error: Could not find \"$_\" column in ${infoFile}!\n" unless exists($colIdx{$_});
    }
    
    
# creating list of exons out of infoFile...
    @{$listref}=();
    %{$queryIdxHash}=() if defined($queryIdxHash);
    # stores the position of query-key in exon-list (optional)
    @{$refexonsList}=() if defined($refexonsList);
    # stores (the position of) reference-exons (IDs) for (in) exon-list (optional)
    @{$exonIDlist}=() if defined($exonIDlist);
    # stores only the Exon-IDs as list (optional)

    my $cPos=0;
    my $colTotal=@tcolumns; #0;
    my $lastRslicePos=undef;
    my $firstRslicePos=undef;
   
    do {
         # while(<INFOF>) {
			if ($headerl =~ m/^[^\#]/) {
  	    chomp($headerl);
	      my @cLi = split(/\s+/, $headerl);
	      #$colTotal = @cLi if(! $cPos); ## first line defines number of columns constraint for rest of lines
	      
			  # add estimated protein length (transcript length divided by 3) if missing protLength column????
		
	      push(@{$listref}, \@cLi ); #print ".";
	     die "ERROR: the number of data columns is not the same as in column-definition (pos${cPos})!\n" unless(@cLi eq $colTotal);
	    
	    $queryIdxHash->{$listref->[$cPos]->[$colIdx{$queryColName}]} = $cPos if defined($queryIdxHash);
	    $firstRslicePos = SynBlast::MyUtils::min($firstRslicePos, $listref->[$cPos]->[$colIdx{'startPos'}]);

	    $lastRslicePos = SynBlast::MyUtils::max($lastRslicePos, $listref->[$cPos]->[$colIdx{'endPos'}]);
	    push(@{$exonIDlist}, $listref->[$cPos]->[$colIdx{$queryColName}]) if defined($exonIDlist);
	    if(defined($refexonsList)) {
		
		if( $listref->[$cPos]->[$colIdx{'geneID'}] eq  $infoH{"focalgene"} ) {
			    # save reference exons corresponding to infofile-ensemblGeneID...
		    push(@{$refexonsList},  $listref->[$cPos]->[$colIdx{$queryColName}]);
		}
	 }
	    ++$cPos;
	}
 } while(defined($headerl=<INFOF>));
    
    close(INFOF);
    my $usedFlankS = SynBlast::MyUtils::max($listref->[$queryIdxHash->{$refexonsList->[0]}]->[$colIdx{'startPos'}] - $listref->[0]->[$colIdx{'startPos'}] + 1,
			       $listref->[-1]->[$colIdx{'endPos'}] - $listref->[$queryIdxHash->{$refexonsList->[-1]}]->[$colIdx{'endPos'}] + 1);
    my $realSlength = $lastRslicePos - $firstRslicePos + 1;
    return(\%colIdx, \%infoH, $usedFlankS, $realSlength);
}



1;
