#!/usr/bin/env perl
# Last changed Time-stamp: <2008-08-19 15:03:56 joe>
# 2006-05-09
# joe@bioinf.uni-leipzig.de
#######################################################

use strict;
use warnings;

package SynBlast::EnsEMBLaccess; 
require Exporter;
my @ISA = qw(Exporter);
my @EXPORT = qw(loadRegistry getDBadaptor storeHomologsForGID checkComparaData getRelMonthString);
my @EXPORT_OK = @EXPORT;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use SynBlast::MyUtils;
use DBI;

############################################################
# queries compara db by reference gene stable id for homologs, and returns success hash 
# for homologs overlapping specified target interval/location
############################################################
sub doCheckComparaByGID {
    my ($geneAdapt, $orgN, $coordS, $seqR, $startP, $endP, $refGeneID, $verboseFl) = @_;
   
    my $refgene = $geneAdapt->fetch_by_gene_stable_id("${refGeneID}");
    die "ERROR: could not fetch_by_gene_stable_id \"${refGeneID}\" (maybe wrong database?)\n" unless(defined $refgene);
    my @refgenes = @{$refgene->get_all_Genes()};
    my %successHash=();
    my $dummy=0;
    foreach my $sgen (@refgenes) {
	next if($sgen->display_id() ne ${refGeneID});
	++$dummy;
	my $homolList = $sgen->get_all_homologous_Genes();
	my $eOrgName= join(" ", split(/_/, $orgN));
	die "ERROR: a second round for same geneID not allowed!\n" unless ($dummy <= 1);
	%successHash=();
	print "############\n"
	    , "Homologous Genes for gene \""
	    , $sgen->display_id()
	    , "\" are:\n" if($verboseFl);
	
	foreach my $hgenRef (@{$homolList}) {
	    my $successFl=0;
	    print $hgenRef->[2] 
		, "\t"
		, $hgenRef->[0]->display_id()
		, "\n" if($verboseFl);
	    if((${eOrgName} eq $hgenRef->[2]) && checkLocation($hgenRef->[0], $orgN, $coordS, $seqR, $startP, $endP, $verboseFl)) {
		++$successFl;
	    }
	    $successHash{$hgenRef->display_id()} = $successFl if($successFl);
	}
	print "#############\n" if($verboseFl); 
    }
    return(\%successHash);
}


sub checkLocation {
    my ($homolGene, $orgN, $coordS, $seqR, $startP, $endP, $verboseFl) = @_;
    return(0) if($homolGene->slice->seq_region_name() ne $seqR);
    return(0) if($homolGene->slice->coord_system_name() ne $coordS);
    my @gtranscripts = @{$homolGene->get_all_Transcripts()};
    my $trans = $gtranscripts[0];
    if($trans->translation()) { #print "coding...\n";
	return(SynBlast::MyUtils::getOverlapSize($homolGene->slice->start() -1 + $trans->coding_region_start(), $homolGene->slice->start() -1 + $trans->coding_region_end(), $startP, $endP));
    }
    return(0);
}

############################################################
# queries compara db by coordinates/positions for homologs, and returns success hash 
# for those genes within interval that are annotated as homologs to the reference gene (ID)
############################################################
sub doCheckComparaByPos {
    my ($slAdapt, $orgN, $coordS, $seqR, $startP, $endP, $refGeneID, $verboseFl) = @_;
    my $slice = undef;
    $slice = $slAdapt->fetch_by_region($coordS, $seqR, $startP, $endP, 1);
    if(defined($slice)) {
	my @sgenes = @{$slice->get_all_Genes()};
	my %successHash=();
	foreach my $sgen (@sgenes) {    
	    my $successFl=0;
	    
	    if($sgen->display_id() eq $refGeneID) { # identity (same gene/organism)
		++$successFl;
		$successHash{$sgen->display_id()} = $successFl if($successFl);
		next;
	    }
	    
	    my $homolList = $sgen->get_all_homologous_Genes();
	    print "############\n"
		, "Homologous Genes for gene \""
		, $sgen->display_id()
		, "\" are:\n" if($verboseFl);
	    foreach my $hgenRef (@{$homolList}) {
		print $hgenRef->[2] 
		    , "\t"
		    , $hgenRef->[0]->display_id()
		    , "\n" if($verboseFl);
		if(defined($refGeneID) && ($hgenRef->[0]->display_id() eq $refGeneID)) {
		    ++$successFl;
		}
	    }
	    print "#############\n" if($verboseFl);
	    $successHash{$sgen->display_id()} = $successFl if($successFl);
	}
	return(\%successHash);
    } else { print "\nCould not retrieve slice!\n"; return(undef); }
}




sub loadRegistry {
    my ($initF, $host, $user) = @_;
    $host = 'ensembldb.ensembl.org' if(! defined($host));
    $user = 'anonymous' if(! defined($user));
    
    if(! defined($initF)) { # no initfile specified -> use automatic retrieval function and ensembl standard web location (or if defined, the selected (local) host instead)
	#print STDERR "Using latest database version from ${host}, user is ${user}.\n";
	Bio::EnsEMBL::Registry->load_registry_from_db( -host => $host,
						       -user => $user,
						       -verbose => "0" );
      } else {
	  Bio::EnsEMBL::Registry->load_all($initF);
	}
}



sub getDBadaptor {
    my ($orgN, $dbType, $adType, $verboseFl) = @_;
    print STDERR ">getting adaptor for (${orgN}, ${dbType}, ${adType})...\n" if($verboseFl);
    my $db_adaptor = Bio::EnsEMBL::Registry->get_adaptor($orgN,$dbType,$adType);
    print STDERR ">APIversion is " . Bio::EnsEMBL::Registry->software_version() . "\n\n" if($verboseFl);
    if(! defined($db_adaptor)) {
	warn "Error: Could not connect to database!\n(Params were: ${orgN}, ${dbType}, ${adType})\n"; 
   }
    return($db_adaptor);
}


############################################################
# updates a hash grouped by targetOrganismName -> refGID -> seqRegName -> values 
# values are containing the members of the homology relations for refGID
############################################################
sub storeHomologsForGID {
    my ($member, $refGID, $targOrg, $hashref, $homologyAdapt) = @_;
    %{$hashref} = () if(! defined($hashref));
    return(-1) if(exists($hashref->{$targOrg}->{$refGID}));
    #%{$hashref->{$targOrg}->{$refGID}} = (); # reset subhash for targOrg->refGID
    my $countr=0;
    # get all the homologies where refGID and targOrg is involved...
	eval {
		$homologyAdapt = SynBlast::EnsEMBLaccess::getDBadaptor('Multi', 'compara', 'Homology') if(! defined($homologyAdapt));
		if(defined($homologyAdapt)) {
			my $homologies=$homologyAdapt->fetch_all_by_Member_paired_species($member, $targOrg);
			print STDERR ">fetching EnsemblCompara-homologs for (GeneID=" . $member->stable_id() . ", targetOrg=${targOrg})...";
			print STDERR ((! @{$homologies}) ? "none" : "ok") . "\n";
			foreach my $homology (@{$homologies}) {
				my $hgenList = $homology->gene_list();
				foreach my $cmember (@{$hgenList}) {
					my $chgenID = $cmember->stable_id();
					my $chgenChr = $cmember->chr_name();
					my $chgenName = $cmember->get_Gene();
					$chgenName = $chgenName->external_name() if(defined($chgenName));
					push(@{$hashref->{$targOrg}->{$refGID}->{$chgenChr}}, [ ($cmember->stable_id()
																		, $cmember->chr_start()
																		, $cmember->chr_end()
																		, $cmember->get_longest_peptide_Member()
																		, $homology->description
																		, $chgenName) ] );
					++$countr;
				}
			} 
		} else { print STDERR "ERROR: fetching homologs for (GeneID=" . $member->stable_id() . ", targetOrg=${targOrg}) failed!!!\n"; } 
	};
	if($@) { print STDERR "ERROR: fetching homologs for (GeneID=" . $member->stable_id() . ", targetOrg=${targOrg}) failed: " . $@ .  "\n"; }
	return($countr);
}


############################################################
# checks the target location on overlapping with (prestored) comparaData, and, if no hit available,
# queries the gene annotations via Ensembl core database for target location
# In both cases, retrieve the protein sequence of target genes, if retrieveFlag set
############################################################
sub checkComparaData {
    my ($comparaData, $slice_adaptor, $orgName, $coordSys, $seqRegName, $fromPos, $toPos, $refGID, $verboseFl, $retrieveFl) = @_;
    # comparaData is specific for orgName and refGID
    my $txtOut="";
    my %hashOut=();
    my %geneNames=();

    my $targetDir = "${refGID}-homologs";
    my $knownOrtDir = "comparaHits";
    my $unknownOrtDir = "unknownHits";

    if($retrieveFl>1) {
	mkdir $targetDir unless(-d $targetDir);
	mkdir "${targetDir}/${knownOrtDir}" unless(-d "${targetDir}/${knownOrtDir}");
	mkdir "${targetDir}/${unknownOrtDir}" unless(-d "${targetDir}/${unknownOrtDir}");
    }
   
   if(! (defined($comparaData) || defined($slice_adaptor))) {
   		### if both first args are undef -> 
   		#### try to fetch all homologs for current refGID now (no prestored comparaData-hash was available)
   		#### and obtain additionally slice_adapt
			SynBlast::EnsEMBLaccess::loadRegistry(); 
			%{$comparaData}=();
		my $member_adaptor = SynBlast::EnsEMBLaccess::getDBadaptor("Multi", "compara", "Member");
		my $homologyAdapt = SynBlast::EnsEMBLaccess::getDBadaptor('Multi', 'compara', 'Homology');
		my $cmember=$member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE', $refGID);
		die "ERROR: could not get member for ${refGID}!\n" unless (defined $cmember);
		  ### query compara db
			SynBlast::EnsEMBLaccess::storeHomologsForGID($cmember, $cmember->stable_id(), $orgName, $comparaData, $homologyAdapt);
			if(defined($comparaData->{$orgName}) && exists($comparaData->{$orgName}->{$refGID})) {
					%{$comparaData}=%{$comparaData->{$orgName}->{$refGID}};
			} else {
					$comparaData=undef;
					$slice_adaptor = SynBlast::EnsEMBLaccess::getDBadaptor($orgName, "core", "Slice");
				}
				
   	}
   

    #print "seqRegName is= ${seqRegName}\n";
    #print "keys of comparaData are=" . join(":", (keys %{$comparaData})) . "\n";
    if(defined($comparaData) && exists($comparaData->{$seqRegName})) {
	#print "## look for overlapp in positions...\n";
	my $successCountr=0;
	my %txtOutLines=();
	foreach my $cHom (@{$comparaData->{$seqRegName}}) {
	    my %curCompH=();
	    my $ptxtOut="";
	    if(SynBlast::MyUtils::getOverlapSize($fromPos, $toPos, $cHom->[1], $cHom->[2]) > 0) {
	      ++$successCountr;
	      $hashOut{$cHom->[4]} = {()} if(! exists($hashOut{$cHom->[4]}));
	#	$txtOut .= "COMPARA-hit:"; #for ${refGID}";
	      $ptxtOut .= $cHom->[4] . ":";
	      my $qGID = $cHom->[0];
	      $geneNames{$qGID}=$cHom->[5] if(defined($cHom->[5]));
	      #$txtResult .= $cHom->[0] . "\n";
	      if($verboseFl) {
		$ptxtOut .= " ${refGID}-ortholog \"${qGID}\" at\t${orgName}\t" 
		    . (defined($coordSys) ? ${coordSys} : "top-level" )
		    . "\t${seqRegName}\t${fromPos}\t${toPos}\n";
	      } else {
		$ptxtOut .= "${qGID}";
	      }
	      
	     if($retrieveFl) {
		
		# save corresponding protein (longest transcript only) of the known homolog targetGeneID
		# if retrieveFl > 1
				
		$ptxtOut .=  "_(";
		my @protList = ();
		my $lPepMember=$cHom->[3];
		my $trans=$lPepMember->get_Transcript;
		my $extname = undef;
		if (my $extn = $trans->external_name()) {
		  $extname = $extn;
		}
		if($trans->translation()) {
		  # test on coding transcript...
		  my $pep=$trans->translate();
			my $curPepID=$pep->display_id();
			my $curPepName=$extname;
			push(@protList, (defined($curPepName) ? "${curPepID} [${curPepName}]" : "${curPepID}"));
			my $hkeyStr="${curPepID}";
			$hkeyStr.= "_${curPepName}" if(defined($curPepName));
			$curCompH{$hkeyStr} = [ ($curPepID, $curPepName) ];
				    
			if($retrieveFl > 1) {
			    my $protFile="${targetDir}/${knownOrtDir}/${orgName}.${curPepID}.fa";
			    if(-f "${protFile}") {
				print STDERR ">\"${protFile}\" did already exist, skip downloading...\n";#if($verboseFl);
			    } else {
				my $curPepLength=$pep->length();
				my $curPepSeq=$pep->seq();
				open(SEQFILE, '>', "${protFile}")  or die "Error while creating ${protFile}: $!\n";
#				print SEQFILE ">${curPepID} ${scoords} ${sseq_region}, " . ($sstart -1 + $trans->start()) . "-" . ($sstart -1 + $trans->end()) . " as translated sequence, length=${curPepLength} strand=" . $trans->strand() . (defined($curPepName) ? " ${curPepName}" : "") . "\n";
				print SEQFILE ">${curPepID} " . "missing coordinates" . " as translated sequence, length=${curPepLength} strand=" . $trans->strand() . (defined($curPepName) ? " ${curPepName}" : "") . "\n";
				print SEQFILE "${curPepSeq}\n";
				close SEQFILE;
			    }
			}
		    }
		$ptxtOut .=  join(",", @protList) . ")\n";
	      } else {
		$ptxtOut .=  "\n";
	      }
	      $hashOut{$cHom->[4]}->{$qGID}= [ values(%curCompH) ];
	      $txtOutLines{$ptxtOut}++;
	   }
	  }
	
	$txtOut .= join("", keys(%txtOutLines));
	return($txtOut, \%hashOut, \%geneNames) if($successCountr);
    } 
    #else 
    ## no comparaHit available -> retrieve overlapping genes 
    #($txtOut, $hashOut) = ("noComparaHit", undef);
    $txtOut .=  "NOMATCH:";
    $txtOut .=  " for ${refGID}" if($verboseFl);
    $hashOut{'noComparaHits'} = {()};
    # should be modified to extract the best matching gene/transcript/protein sequence with respect to the blast-hits (putative exon regions) in this region
    ####

    if($retrieveFl) {
	if(defined($slice_adaptor)) {
		my $slice = undef;
		
		#print "fetching by region params=" .  join(" " , ($coordSys, $seqRegName, $fromPos, $toPos)) . "\n";
		$slice = $slice_adaptor->fetch_by_region($coordSys, $seqRegName, $fromPos, $toPos, 1);
		
		if(defined($slice)) {
			my $scoords = $slice->coord_system()->name();
			my $sseq_region = $slice->seq_region_name();
			my $sstart = $slice->start();
			my $send = $slice->end();
			#my $sstrand = $slice->strand();
			my @sgenes = @{$slice->get_all_Genes()};
			$txtOut .=  "\n#There are " . (@sgenes) . " overlapping genes:\n" if($verboseFl);
			my @curNoCompH=();
			foreach my $sgen (@sgenes) {
			#      $txtOut .=  "# " . $sgen->display_id() . " (";
				$txtOut .=  $sgen->display_id() . "(";
				$geneNames{$sgen->display_id()}=$sgen->external_name() if(defined($sgen->external_name()));
		my @protList = ();
	      my @gtranscripts = @{$sgen->get_all_Transcripts()};
#    my $transcount = @gtranscripts;
	      foreach my $trans (@gtranscripts) {
		my $extname = undef;
		    if (my $extn = $trans->external_name()) {
			$extname = $extn;
		    }
		    if($trans->translation()) {
			# test on coding transcript...
			#	if(($trans->coding_region_start() > 0) && ($trans->coding_region_end() <= $slice->length())) { # coding region is completely  within slice???
			my $pep=$trans->translate();
			my $curPepLength=$pep->length();
			my $curPepID=$pep->display_id();
			my $curPepName=$extname;
			my $curPepSeq=$pep->seq();
			my $protFile="${targetDir}/${unknownOrtDir}/${orgName}.${curPepID}.fa";
			push(@protList, (defined($curPepName) ? "${curPepID} [${curPepName}]" : "${curPepID}")); 
			push(@curNoCompH, [ ($curPepID, $curPepName)]);
			if($retrieveFl > 1) {
			    if(-f "${protFile}") {
				print STDERR ">\"${protFile}\" did already exist, skip downloading...\n"; # if($verboseFl);
			    } else {
				open(SEQFILE, '>', "${protFile}")  or die "Error while creating ${protFile}: $!\n";
				print SEQFILE ">${curPepID} ${scoords} ${sseq_region}, " . ($sstart -1 + $trans->start()) . "-" . ($sstart -1 + $trans->end()) . " as translated sequence, length=${curPepLength} strand=" . $trans->strand() . (defined($curPepName) ? " ${curPepName}" : "") . "\n";
				print SEQFILE "${curPepSeq}\n";
				close SEQFILE;
			    }
			    #}
			}
		    }
		}
		$txtOut .=  join(",", @protList) . ");";
		$hashOut{'noComparaHits'}->{$sgen->display_id()}=\@curNoCompH;
	    }
	} else {
	    print STDERR "\n>Could not get slice! (" .  join(", "  , ($coordSys, $seqRegName, $fromPos, $toPos)) . ")\n";
	}
      }
    }
    else {
      $txtOut .=  "\n";
      
    }
    
    return($txtOut, \%hashOut, \%geneNames);
  }






### to automatically retrieve the ensembl website archive version string (monthYear) for a given release number
sub getRelMonthString {
	 my ($queryRelV, $host, $port, $user, $pass, $verbose, $db_version) =
			SynBlast::MyUtils::rearrange( [qw(FOR_RELEASE HOST PORT USER PASS VERBOSE DB_VERSION )], @_);
	
	$queryRelV= Bio::EnsEMBL::Registry->software_version() unless($queryRelV);
	$db_version= Bio::EnsEMBL::Registry->software_version() unless($db_version);
	$host ="ensembldb.ensembl.org" unless($host);
	$user= "anonymous" unless($user);
	$pass="" unless(defined($pass));
	
	if(!defined($port)){
	$port   = 3306;
		if($host eq "ensembldb.ensembl.org"){
			if( !defined($db_version) or $db_version >= 48){
				$port = 5306;
			}	
		}
	}

	my $db = DBI->connect( "DBI:mysql:host=$host;port=$port" , $user, $pass );
        if(($db_version == 47)  || ($db_version == 48)) {
		$db_version = 49;  ## workaround for missing 'ensembl_website_47' and '..._48' databases...
	}
	my $res = $db->selectall_hashref(   #"use ensembl_website_${db_version}; 
					"select release_id,archive from ensembl_website_${db_version}.ens_release", "release_id");
	
	return($res->{$queryRelV}->{'archive'});
	
}



############################################################


__END__

=head1 SynBlast::EnsEMBLaccess - Wrapper modules for EnsemblAPI specific/related functions



=head1 SYNOPSIS

 require Bio::EnsEMBL::DBSQL::DBAdaptor;
 require Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
 require SynBlast::MyUtils;
 
=head1 DESCRIPTION

functions providing easy access to Ensembl databases

=head1 AUTHOR

Joerg Lehmann, <lt>Ejoe@bioinf.uni-leipzig.deE<gt>

=cut
