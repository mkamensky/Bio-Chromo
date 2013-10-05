package Bio::Chromo::CDS;
# ABSTRACT: Class representing one C<CDS>

use 5.10.0;
use integer;
use overload '""' => '_as_string';
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

=func DataVersion

Version of the data structure, to use when storing. Update when changing the 
underlying data structure.

=cut

sub DataVersion { 1.3 }

=method new

    $cds = new Bio::Chromo::CDS $chromo, $feat;

Create a new L<Bio::Chromo::CDS> object, belonging to chromosome $chromo, and 
representing feature $feat. $feat should be something that looks like a 
L<Bio::SeqFeatureI>.

=cut

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

=method encode

    $str = $self->encode($seq)

This method encodes a string of bases into a shorter string, where each codon 
is represented by one character. Thus, each character in the result 
determines an amino acid. This is a class method.

=cut

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

=method decode

    $seq = $self->decode([$start[,$len[,$str]]]);

Decode a string as produced by L</encode> back into a sequence of bases. $str 
is the string to decode. If not given, the internal string stored in $self is 
used (if it is given, this can be used as a class method). $start and $len 
give the starting index and length of the fragment to decode. They default to 
0 and 1, respectively.

=cut

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

sub _as_string { $_[0]->{'gene'} }

sub _loc {
    defined $_[1] ? $_[0]->{'loc'}[$_[1]] : $_[0]->{'loc'}
}

=method find_exon

    $exon = $cds->find_exon($i);

Find to which exon a given chromosome position belongs. The result is an 
index into the list of exons (locations) for this gene.

=cut

sub find_exon {
    my ($self, $i) = @_;
    for my $n ( 0..$#{$self->_loc} ) {
        return $n if $self->_loc($n)->contains($i);
    }
    ()
}

=method gene_pos

    $pos = $cds->gene_pos($i);

Given an absolute position $i into the chromosome, returns the same position 
in the coordinates of the gene to which it belongs. Honours direction. If $i 
is not within the gene, an empty list is returned.

=cut

sub gene_pos {
    my ($self, $i) = @_;
    my $e = $self->find_exon($i);
    return () unless defined $e;
    my $p = $self->{'lens'}[$e];
    $p + ($self->{'reverse'} ? $self->_loc($e)->end - $i 
                             : $i - $self->_loc($e)->start)
}

=method codon

     $codon = $cds->codon($i);

Return the codon at position $i. Returns a string of 3 base letters. The 
position $i is an absolute position in the chromosome.

=cut

sub codon {
    my ($self, $i) = @_;
    map { $self->decode($_/3) } $self->gene_pos($i)
}

sub _trans {
    require Bio::Perl;
    Bio::Perl::translate_as_string($_[1])
}

=method amino

    $amino = $self->amino($i);

Returns the amino acid coded by codon at location $i. The location is 
interpreted as in L</codon>.

=cut

sub amino {
    my ($self, $i) = @_;
    map { $self->_trans($_) } $self->codon($i)
}

=head1 SEE ALSO

L<Bio::Chromo>

=cut

1;

