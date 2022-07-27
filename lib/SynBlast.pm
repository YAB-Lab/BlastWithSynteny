package SynBlast;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our @EXPORT_OK = qw($VERSION $url);
our @EXPORT = @EXPORT_OK;

our $VERSION = '0.03';
our $url = "http://www.bioinf.uni-leipzig.de/Software/SynBlast/";

# Preloaded methods go here.

1;
__END__

=head1 NAME

SynBlast - a Perl package for Assisting the Analysis of Conserved Synteny Information

Copyright (C) 2008  Joerg Lehmann

=head1 SYNOPSIS

  use SynBlast;

=head1 DESCRIPTION


=head2 EXPORT

None by default.

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
L<getEnsemblProteins.pl>
L<doBlastJobs.pl>
L<doSyntenyFilter.pl>,
L<doSyntenyAlignment.pl>,
L<getBestBlastHit.pl>

=cut
