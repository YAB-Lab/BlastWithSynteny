#!/usr/bin/env perl
# Last changed Time-stamp: <2008-08-19 15:04:13 joe>
# joe@bioinf.uni-leipzig.de
########################################################################

use strict;
use warnings;

package SynBlast::MyUtils;
require Exporter;
my @ISA = qw(Exporter);
my @EXPORT = qw(getSgn formatMyNumber comparaRes2HTML printHash printHashHash printListHash printListList printHashListList printHashListHash min max swapArray printHashList printHashHashList printListHashHash printHashHashListHash printHashHashListList getLeadingZeroStr getOverlapSize swapScalar convertListList checkStrEnd checkRelDir);
 my @EXPORT_OK = @EXPORT; 



=head2 rearrange

 Usage     : rearrange( array_ref, list_of_arguments)
 Purpose   : Rearranges named parameters to requested order.
 Example   : use MyUtils qw(rearrange);
           : rearrange([qw(SEQUENCE ID DESC)],@param);
           : Where @param = (-sequence => $s, 
	         :                 -id       => $i, 
	         :	               -desc     => $d);
 Returns   : @params - an array of parameters in the requested order.
           : The above example would return ($s, $i, $d)
 Argument  : $order : a reference to an array which describes the desired
           :          order of the named parameters.
           : @param : an array of parameters, either as a list (in
           :          which case the function simply returns the list),
           :          or as an associative array with hyphenated tags
           :          (in which case the function sorts the values 
           :          according to @{$order} and returns that new array.)
	   :	      The tags can be upper, lower, or mixed case
           :          but they must start with a hyphen (at least the
           :          first one should be hyphenated.)
 Source    : This function was taken from CGI.pm, written by Dr. Lincoln
           : Stein, and adapted for use in Bio::Seq by Richard Resnick and
           : then adapted for use in Bio::Root::Object.pm by Steve A. Chervitz.
           : This has since been adapted as an exported static method in this 
             class Bio::EnsEMBL::Utils::Argument 
 Comments  : (SAC)
           : This method may not be appropriate for method calls that are
           : within in an inner loop if efficiency is a concern.
           :
           : Parameters can be specified using any of these formats:
           :  @param = (-name=>'me', -color=>'blue');
           :  @param = (-NAME=>'me', -COLOR=>'blue');
           :  @param = (-Name=>'me', -Color=>'blue');
           : A leading hyphenated argument is used by this function to 
           : indicate that named parameters are being used.
           : Therefore, a ('me', 'blue') list will be returned as-is.
           :
           : Note that Perl will confuse unquoted, hyphenated tags as 
           : function calls if there is a function of the same name 
           : in the current namespace:
           :    -name => 'foo' is interpreted as -&name => 'foo'
           :
           : For ultimate safety, put single quotes around the tag:
	         :    ('-name'=>'me', '-color' =>'blue');
           : This can be a bit cumbersome and I find not as readable
           : as using all uppercase, which is also fairly safe:
           :    (-NAME=>'me', -COLOR =>'blue');
	         :
           : Personal note (SAC): I have found all uppercase tags to
           : be more managable: it involves less single-quoting,
           : the code is more readable, and there are no method naming 
           : conlicts.
           : Regardless of the style, it greatly helps to line
	         : the parameters up vertically for long/complex lists.

=cut


sub rearrange {
  my $order = shift;
  $order = shift if($order eq "SynBlast::MyUtils"); #skip object if one provided

  # If we've got parameters, we need to check to see whether
  # they are named or simply listed. If they are listed, we
  # can just return them.
  return @_ unless (@_ && $_[0] && substr($_[0], 0,1) eq '-');

  # Convert all of the parameter names to uppercase, and create a
  # hash with parameter names as keys, and parameter values as values
  my $i = 0;
  my (%param) = map {if($i) { $i--; $_; } else { $i++; uc($_); }} @_;

  # What we intend to do is loop through the @{$order} variable,
  # and for each value, we use that as a key into our associative
  # array, pushing the value at that key onto our return array.
  return map {$param{uc("-$_")}} @$order;
}




sub getSgn {
    my($inp)=@_;
    return(-1) if($inp < 0);
    return(1);
}


sub formatMyNumber {
	my($number, $sepCh, $sepL) = @_;
	return("") if(! defined($number));
	my $numbers = "" . int($number);
	if($sepL > 0) {
	    my $result="";
	    my @tmpA=();
			my $li=length(${numbers});
			while($li> $sepL) {
				push(@tmpA, substr($numbers, $li - $sepL, $sepL));
				$li -= $sepL;
			} 		
			$result = substr($numbers, 0, $li);
		  $result .= ($sepCh . join($sepCh, reverse @tmpA)) if(@tmpA);
		  return($result);
  } else {
  	return($numbers);
  }
}

sub comparaRes2HTML {
    my ($resHash, $compOrg, $ensemblSite, $gID2Name, $protFl) = @_;
    
    %{$gID2Name}=() if(! defined($gID2Name));
    $protFl=0 if(! defined($protFl));
    my @resType= keys %{$resHash};
    return("noResult") unless(@resType);
    warn "Note: comparaResult did not give exactly one result type (" . join(" and ", @resType) . " for target ${compOrg})!\n" unless(@resType == 1);
    my %ensemblRefTypes = ( 'protID' => 'protview?peptide=',
			    'geneID' => 'geneview?gene=',
			    'transID' => 'transview?transcript=');
		
    $ensemblSite = "http://www.ensembl.org" if(! defined($ensemblSite));

    my $compResHtml="";
    
    for(my $resTidx=0; $resTidx < @resType; ++$resTidx) {
	    
	$compResHtml.= $resType[$resTidx] . ": ";
	my @tmpGn=keys %{$resHash->{$resType[$resTidx]}};
	for(my $gidx=0; $gidx < @tmpGn; ++$gidx) {
		$compResHtml .= "<a href=\"" 
			. "${ensemblSite}/${compOrg}/" 
			. ($ensemblRefTypes{'geneID'})
			. $tmpGn[$gidx] 
			.  "\" target=\"_ensembl\">"
			. ( (exists($gID2Name->{$tmpGn[$gidx]})) ? $gID2Name->{$tmpGn[$gidx]} :  $tmpGn[$gidx] )
			. "</a>\n";
		my @tmpPn=@{$resHash->{$resType[$resTidx]}->{$tmpGn[$gidx]}};
		if(@tmpPn && $protFl) {
			$compResHtml .= " (";
			for(my $pidx=0; $pidx < @tmpPn; ++$pidx) {
				my $lire= $tmpPn[$pidx];
				$compResHtml .= "<a href=\"" 
					. "${ensemblSite}/${compOrg}/" 
					. ($ensemblRefTypes{'protID'})
					. $lire->[0]
					.  "\" target=\"_ensembl\">"
					. $lire->[0];
				$compResHtml .= " [" . $lire->[1] . "]" if($lire->[1]);
				$compResHtml .= "</a>\n";
				$compResHtml .= ", " if($pidx + 1 < @tmpPn);
			}
			$compResHtml .= ")";
		}
		$compResHtml .= ";<br/> " if($gidx + 1 < @tmpGn);
	}
	$compResHtml .= "<br/>" if($resTidx + 1 < @resType);
     }
     $compResHtml .= "\n";	
     return($compResHtml);
}
    


sub printHash {
    my ($hashref, $sep1, $sep2, $preText) = @_;
    if (! defined($preText)) {
	$preText=""; }
    else {
	$preText=$preText . $sep1;
    }
    foreach my $curE (sort keys %$hashref) {
	print $preText . $curE . $sep1 . $hashref->{$curE} . $sep2;
    }
}



sub printHashHash {
    my ($hashref, $sep1, $sep2) = @_;
    foreach my $chr (sort keys %$hashref) {
	printHash(\%{$hashref->{$chr}}, $sep1, $sep2, $chr);
	print $sep2;
    }
}


sub printListHash {
    my ($listref, $sep1, $sep2, $sep3, $preText) = @_;
    if (! defined($preText)) {
	$preText=""; 
    } else {
	$preText=$preText . $sep1;
    }
    for(my $chri=0; $chri < @$listref; $chri++) {
	printHash(\%{$listref->[$chri]}, $sep1, $sep2, $preText . $chri);
    }
}


sub printHashList {
    my ($hashref, $sep1, $sep2) = @_;
    foreach my $curE (keys %$hashref) {
	print $curE . $sep1 . join($sep1, @{$hashref->{$curE}} ) . $sep2;
    }
}


sub printListList {
    my ($lref, $sep1, $sep2, $preText) = @_;
    if(! defined($lref)) {
	return(0);
    }
    if (! defined($preText)) {
	$preText=""; 
    } else {
	$preText=$preText . $sep1;
    }
    foreach my $curEE (@{$lref}) {
	if(! defined($curEE)) {
	    return(0);
	}
	print $preText . join($sep1, @$curEE) . $sep2;
    }
    return(1);
}


sub printHashHashList {
    my ($hashref, $sep1, $sep2, $sep3) = @_;
    foreach my $chr (sort keys %$hashref) {
	print "key $chr: " . (keys %{$hashref->{$chr}}) . $sep3;
	printHashList(\%{$hashref->{$chr}}, $sep1, $sep2);
    }
}


sub printListHashList {
    my ($listref, $sep1, $sep2, $sep3) = @_;
    for(my $li=0; $li<@{$listref}; $li++) {
	print "idx ${li}: " . (keys %{$listref->[$li]}) . $sep3;
	printHashList(\%{$listref->[$li]}, $sep1, $sep2);
    }
}



sub printListHashHash {
    my ($listref, $sep1, $sep2, $sep3) = @_;
    for(my $li=0; $li<@{$listref}; $li++) {
	print "idx ${li}: " . (keys %{$listref->[$li]}) . $sep3;
	printHashHash(\%{$listref->[$li]}, $sep1, $sep2);
	print $sep3;
    }
}



sub printHashListList {
    my ($hashref, $sep1, $sep2, $sep3) = @_;
    foreach my $curE (sort keys %$hashref) {
	print "key $curE: " . @{$hashref->{$curE}} . $sep3;
	printListList(\@{$hashref->{$curE}}, $sep1, $sep2);
    }
}


sub printHashHashListList {
    my ($hashref, $sep1, $sep2, $sep3, $sep4) = @_;
    foreach my $curE (sort keys %{$hashref}) {
	print "key $curE: " . (keys %{$hashref->{$curE}}) . $sep4;
	printHashListList(\%{$hashref->{$curE}}, $sep1, $sep2, $sep3);
    }
}


sub printHashHashListHash {
    my ($hashref, $sep1, $sep2, $sep3, $sep4) = @_;
    foreach my $curE (sort keys %{$hashref}) {
	print "key $curE: " . (keys %{$hashref->{$curE}}) . $sep4;
	printHashListHash(\%{$hashref->{$curE}}, $sep1, $sep2, $sep3);
    }
}





sub printHashListHash {
    my ($hashref, $sep1, $sep2, $sep3) = @_;
    if(! defined(keys %$hashref)) {
	return(0);
    }
    foreach my $curE (sort keys %$hashref) {
	print "key $curE: " . @{$hashref->{$curE}} . $sep3;
	printListHash(\@{$hashref->{$curE}}, $sep1, $sep2, "", $curE);
    }
}



sub min { #Numbers.
	  my $min = shift;
    while((! defined($min)) && (@_ > 0)) { $min = shift; }
    foreach ( @_ ) { $min = $_ if(defined($_) && ($_ < $min)); }
    return $min;
}


sub max { #Numbers.
    my $max = shift;
    while((! defined($max)) && (@_ > 0)) { $max = shift; }
    foreach ( @_ ) { $max = $_ if(defined($_) && ($_ > $max));	}
    return $max;
}



sub swapArray {
    my ($aref, $idx1, $idx2) = @_;
    my $tmp = $aref->[$idx1];
    $aref->[$idx1]=$aref->[$idx2];
    $aref->[$idx2]=$tmp;
}


sub getLeadingZeroStr {
    my ($number, $totalL) = @_;
    return(("0" x (${totalL}-length("${number}"))) . "${number}");
}


sub getOverlapSize {
    my ($s1start, $s1end, $s2start, $s2end) = @_;
    ($s1start <= $s1end) or swapScalar(\$s1start, \$s1end);
  
    ($s2start <= $s2end) or swapScalar(\$s2start, \$s2end);

    return(0) if($s1end < $s2start);
    return(0) if($s2end < $s1start);
    return(  min($s1end, $s2end) - max($s2start, $s1start) + 1);
}  



sub swapScalar {
    my ($ref1, $ref2) =@_;
    my $tmpref=$ref2;
    $ref2=$ref1;
    $ref1=$tmpref;
}



sub convertListList {
    my($llistref, $outref) = @_;
    @{$outref}=();
    for(my $lid=0; $lid<@{$llistref}; $lid++) {
	for(my $ld=0; $ld< @{$llistref->[$lid]}; $ld++) {
	    $outref->[$ld]->[$lid]= $llistref->[$lid]->[$ld];
	}
    }
    return($outref);
}



sub checkStrEnd { # removes last char of string, if ending with $char
    my($strref, $char) = @_;
    die "ERROR: char needed as 2nd argument (checkStrEnd)\n" unless (length($char) == 1);
    my $actionC=0;
    if(substr(${$strref}, -1) eq "${char}" ) {
	${$strref} = substr(${$strref}, 0, (length(${$strref}) -1));
	$actionC++;
    }
    return($actionC);
}



sub checkRelDir { #append dirstring to basedirString, if relative (not starting with "/")
     my($strref, $basedirStr) = @_;
     if(substr(${$strref}, 0, 1) ne "/") {
	 ${$strref} = "${basedirStr}/${$strref}";
	 return 1;
     }
     return 0; #str was already non-relative
 }


1;
