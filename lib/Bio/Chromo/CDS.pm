package Bio::Chromo::CDS;
# ABSTRACT: Class representing one C<CDS>

use 5.10.0;
use integer;
use overload '""' => 'as_string';
use warnings;
use strict;
# This is required here, even though we don't create them explicitly, since 
# we load them from the yaml file
use Bio::Location::Simple;

=head1 SYNOPSIS

    # get a Bio::Chromo object $chromo from somewhere
    my $cds = $chromo->[0]; # first CDS

    # $i is a nucleotide position in the chromosome

    # which exon?
    my $ex = $cds->find_exon($i);
    # what is the position in the gene?
    my $pos = $cds->gene_pos($i);
    # what codon are we in?
    my $codon = $cds->codon($i);
    # what amino acid we code?
    my $amino = $cds->amino($i);

=cut

sub DataVersion { 1.2 }

sub new {
    my ($class, $chromo, $f) = @_;
    my @loc = $f->location->each_Location;
    unless ( @loc ) {
        $chromo->_warn(
            "Feature at " . $f->start .  " has no location, ignoring"
        );
        return ()
    }
    my %r = (
        start => $f->start,
        end => $f->end,
        reverse => ($f->strand == -1),
        loc => \@loc,
    );
    my @g = $f->get_tag_values('gene');
    $chromo->_warn("More than one gene: @g") unless @g == 1;
    $r{'gene'} = shift @g;
    # store for each exon the total length of the gene up to it, according to 
    # the strand
    my @l = (0);
    my $seq = '';
    # unfortunately, Bio::PrimarySeq::subseq() is highly inefficient...
    my $fullseq = $f->entire_seq->{'seq'};
    for my $ex ( $r{'reverse'} ? reverse @loc : @loc ) {
        push @l, $l[-1] + $ex->length;
        #my $sseq = $f->entire_seq->subseq($ex->start, $ex->end);
        my $sseq = substr($fullseq, $ex->start-1, $ex->end-$ex->start+1);
        $seq .= $r{'reverse'} ? reverse $sseq : $sseq;
    }
    pop @l;
    if ( $r{'reverse'} ) {
        $r{'lens'} = [reverse @l];
        my $n = $seq =~ tr/ACGT/TGCA/;
        $chromo->_warn("Incorrect number of bases: $n vs. "  . length($seq)) 
            unless $n==length($seq);
    } else {
        $r{'lens'} = \@l;
    }
    my $self = bless \%r => $class;
    $r{'seq'} = $self->encode($seq);
    $self
}

{

my @num2base = qw(A C G T);
my %base2num = map { $num2base[$_] => $_ } 0..$#num2base;

sub encode {
    my ($self, $seq) = @_;
    my $vec = '';
    my $l = length($seq);
    my $i = $l/3*4;
    while ( $l > 2 ) {
        vec($vec, --$i, 2) = 1;
        vec($vec, --$i, 2) = $base2num{substr($seq, --$l, 1)};
        vec($vec, --$i, 2) = $base2num{substr($seq, --$l, 1)};
        vec($vec, --$i, 2) = $base2num{substr($seq, --$l, 1)};
    }
    $vec
}

sub decode {
    my ($self, $start, $len, $vec) = @_;
    $vec //= $self->{'seq'};
    $start //= 0;
    $len //= 1;
    my @res;
    for my $i ( $start..$start+$len-1 ) {
        push @res, join('', map { $num2base[vec($vec, $_, 2)] } 4*$i..4*$i+2);
    }
    @res
}
}

sub as_string { $_[0]->{'gene'} }

sub loc {
    defined $_[1] ? $_[0]->{'loc'}[$_[1]] : $_[0]->{'loc'}
}

sub find_exon {
    my ($self, $i) = @_;
    for my $n ( 0..$#{$self->loc} ) {
        return $n if $self->loc($n)->contains($i);
    }
    ()
}

# return the position of chromosome position $i inside this gene. Honours
# direction.
sub gene_pos {
    my ($self, $i) = @_;
    my $e = $self->find_exon($i);
    return () unless defined $e;
    my $p = $self->{'lens'}[$e];
    $p + ($self->{'reverse'} ? $self->loc($e)->end - $i 
                             : $i - $self->loc($e)->start)
}

sub codon {
    my ($self, $i) = @_;
    map { $self->decode($_/3) } $self->gene_pos($i)
}

sub trans {
    require Bio::Perl;
    Bio::Perl::translate_as_string($_[1])
}

sub amino {
    my ($self, $i) = @_;
    map { $self->trans($_) } $self->codon($i)
}

=head1 SEE ALSO

L<Bio::Chromo>

=cut

1;

