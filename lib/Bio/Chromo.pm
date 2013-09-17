package Bio::Chromo;
# ABSTRACT: List of genes in a chromosome

use 5.10.0;
use strict;
use warnings FATAL => 'all';

use base qw(List::Sorted);
use File::Spec::Functions qw(catfile catdir splitpath);

=head1 SYNOPSIS

    $chromo = new Bio::Chromo 'NC_000002';
    @aminos = $chromo->aminos(272244); # D, T
    @codons = $chromo->codons(272244); # gac, acc

=head1 DESCRIPTION

This class represents the list of C<CDS> features in a chromosome.  
The list is held sorted (by inheriting from L<List::Sorted>) for speed of 
lookup. The main purpose of this class is currently to easily retrieve the 
amino acid produced by a particular location in the chromosome. This is 
achieved through the L</aminos()> method. The other notable feature of this 
classes is that it automatically caches the chromosome data it needs, which 
drastically reduces the time, memory and bandwidth consumed.

=head1 BIOLOGICAL BACKGROUND

This section contains a short (and imprecise) summary of the biological 
background required to understand the problem solved by this module.

The DNA consists of a fixed (per species) number of I<chromosomes>. Each 
chromosome can be dealt with independently, and so we always consider a fixed 
chromosome. In humans, the chromosomes are enumerated 1-22, X, Y. The 
chromosome is a sequence of the letters A,T,C,G. Each such letter is called a 
I<base> or a I<nucleotide>. The letters A and T and the letters C and G are 
dual to each other. The process of turning each letter to its dual in a given 
sequence is termed I<complementing> the sequence.

The purpose of the DNA is to code (and produce) I<proteins>. A protein is 
coded by a segment of the chromosome called a I<gene>. The protein itself is 
a sequence of I<amino acids>, and the protein is coded by coding this 
sequence. Each kind of amino acid is also symbolised by a fixed latin letter.  
An amino acid is coded by a sequence of 3 bases in the gene. Such a sequence 
is called a I<codon>. Thus, the main goal of this module is to compute in 
which codon a particular base lies. We further have the following 
terminology:

=over

=item Feature

Generic name for an interesting part of the chromosome. Types of features are 
identified by their I<primary tag>. For us, it seems that the only one of 
interest is B<CDS>.

=item CDS

Essentially, the same as a gene. However, every gene may have more than one 
interpretation, in the sense that the same (approximate) region in the 
chromosome might code several different proteins. These different variants 
correspond to different CDS features, though they are considered to be the 
same gene. It is impossible to determine, from the chromosome data alone, 
which interpretation is actually used: This depends on biological input 
external to the DNA (and might change, e.g., from one cell to another).  
Therefore, we consider all of them.

In this module, a CDS feature is represented by a helper class, 
L<Bio::Chromo::CDS>. Each L<Bio::Chromo> object is essentially a sorted array 
of L<Bio::Chromo::CDS> objects.

=item Gene

Largest meaningful unit. Each gene codes a protein, though there might be 
several variants, as explained above. Each gene has a name, like 'FAM110C', 
and these names are different for different genes on the same chromosome 
(though not for different CDSs that correspond to the same gene). 

The part of the gene that actually corresponds to the protein need not be 
contiguous: each contiguous part is called an I<exon>. The part that actually 
codes the protein is obtained by concatenating the exons in one gene, and 
then possibly reversing (the order) and complementing the resulting sequence.

The list and boundaries of the exons might differ in different CDSs 
corresponding to the same gene. For this reason, we are not directly 
interested in the genes, and they have no counterpart in the code.

=item Exon

As explained above, this is a connected component of the part of a gene that 
actually codes a protein. Each of this is represented in the code by a 
L<Bio::Location::Simple> object.

=back

=head2 Algorithm for finding the codon

Summarising all the information above, we discover that to find the codons of 
the base at position N, is logically equivalent to:

=over

=item 1.

Find all CDS features containing position N

=item 2.

For each such feature, concatenate the exons to obtain a string S

=item 3.

Compute the location M of position N within the string S

=item 4.

If the direction (strand) of the gene is reversed, reverse and complement S, 
updating M.

=item 5.

Compute L=M/3, rounding down. The codon consists of the three bases within S 
starting at 3*L.

=back

In reality, we perform a slightly different computation, and also compute 
directly the amino acid, rather than the codon (the translation from codons 
to amino acids is completely determined). This is done for reasons of 
efficiency, and should be equivalent.

=cut

=method new()

    my $chromo = new Bio::Chromo [B<$seq>]

Create a new L<Bio::Chromo> object. The B<$seq> argument can be either a 
L<Bio::SeqI> type object (or anything else that has a I<get_SeqFeatures> 
method), or an id, such as C<NC_000002>. The object is initialised using 
either L<init_from_seq()> or L<init_from_id()>, as appropriate.

=cut

sub new {
    my ($class, $seq) = @_;
    my $self = $class->SUPER::new;

    if ( $seq ) {
        if ( eval { $seq->can('get_SeqFeatures') }) {
            $self->init_from_seq($seq);
        } else {
            $self->init_from_id($seq);
        }
    }
    $self
}

# order the CDS first by the starting point, and within each variation by the
# end point.
sub cmp {
    my $self = shift;
    $_[0]->{'start'} <=> $_[1]->{'start'} or $_[0]->{'end'} <=> $_[1]->{'end'}
}

=method aminos()

    @aminos = $self->aminos($i)

Return the list of possible amino acids coded by a possible codon to which 
the given base $i belongs. Thus, each element in the result is one upper case 
letter, or C<*> for the end of the sequence.

=method codons()

    @codons = $self->codons($i)

Return the list of possible codons to which  the given base $i belongs. Thus, 
each element in the result consists of a string of length 3 on the alphabet 
C<ACGT>.

=cut

for my $met ( qw(amino codon) ) {
no strict 'refs';
my $mets = $met . 's';
*$mets = sub {
    my ($self, $i) = @_;
    map { $self->[$_]->$met($i) } $self->genes($i)
}
}

=method genes()

    @genes = $self->genes($i);

Returns the indices of all genes containing the given position. More 
precisely, each element returned is an index into $self of a feature (C<CDS>) 
containing absolute position $i on the sequence. Note that the direction 
(strand) plays no role here. The features themselves, as L<Bio::Chromo::CDS> 
objects, can be accessed by usual array notation into $self, as in

    for my $g ( @genes ) {
        my $feat = $self->[$g];
        # $feat is a Bio::Chromo::CDS object now
        # ...
    }

=cut

sub genes {
    my ($self, $i) = @_;
    $self->_croak("Called 'genes' without an index") unless defined $i;
    my $f = { start => $i+1, end => $i };
    my $p = $self->find_pos($f);
    my @res;
    while ( --$p >= 0 ) {
        my $c = $self->[$p];
        last if $i > $c->{'end'};
        push @res, $p;
    }
    @res
}

=method init_from_seq()

    $self->init_from_seq($seq)

Initialise the object from the given L<Bio::SeqI> object C<$seq>. In fact, 
C<$seq> is only required to have a I<get_SeqFeatures> method, which returns a 
list of objects that look like L<Bio::SeqFeatureI>. The resulting object will 
hold a list of items corresponding to the B<CDS> items in C<$seq>.

Returns true iff the initialisation was successful

=cut

sub init_from_seq {
    my ($self, $seq) = @_;
    for my $feat ($seq->get_SeqFeatures('CDS')) {
        my $cd = $self->_FeatClass->new($self, $feat);
        next unless defined $cd;
        #$self->_info("  CDS $cd");
        $self->insert($cd);
    }
    1
}

=method init_from_id()

    $self->init_from_id($id)

Initialise the object from the sequence with id C<$id>, a string that looks 
something like B<NC_000002>. This is logically equivalent to fetching the 
sequence given by C<$id> from GenBank, and then calling L<init_from_seq()> 
with it, except we store a cache of the actual structures, and so this is 
potentially much faster and less memory and network consuming. Both the cache 
file and the GenBank file will be stored on the local machine, in the 
directory given by L<DataDir()>. These files need to be erased if the 
GenBank file should be re-downloaded from the site (for example, if there is 
a new version on the web).

Returns true iff the initialisation was successful

=cut

sub init_from_id {
    my ($self, $id) = @_;
    my $dir = $self->DataDir;
    my $dver = $self->DataVersion;
    my $cache = catfile($dir, "$id.yml");
    my $gb = catfile($dir, "$id.gbk");
    my $seq;
    if ( -r $gb ) {
        # We already have a genbank file
        $self->_info("Trying to load cache from '$cache'");
        if ( my $err = $self->load($cache, $dver) ) {
            # gb exists, but not cache. Get a seq
            $self->_warn("$err\n  ...(Re)generating from '$gb'");
            require Bio::SeqIO;
            my $stream = new Bio::SeqIO -file => $gb;
            $self->_croak("Failed to create stream") unless $stream;
            $seq = $stream->next_seq;
            $self->_info("Loaded sequence from $gb");
        } else {
            return 1;
        }
    } else {
        # we need to fetch the sequence from the web
        $self->_info(
            "Genbank file '$gb' does not exist, trying to retrieve...");
        $seq = $self->_fetch_gb($id, $gb, $dir);
    }
    return () unless defined $seq;
    my @ret = $self->init_from_seq($seq);
    # write yaml cache
    $self->_info("  Writing cache to '$cache'");
    if ( my $err = $self->store($cache, $dver) ) {
        $self->_warn($err);
    }
    @ret
}

=method DataDir()

Returns the directory where the genbank files and the cache files are stored.  
Taken from the environment variable C<$GENBANK_DIR>, or C<$HOME/genbank> by 
default.

=cut

our $DataDir = $ENV{'GENBANK_DIR'} || catdir($ENV{'HOME'}, 'genbank');

sub DataDir { $DataDir }

=method DataVersion()

Version of the data stored in the cache. Should be increased whenever the 
structure of the L<Bio::Chromo::CDS> package is changed, or the format of the data 
stored in the C<yaml> cache files is otherwise modified.

=cut

sub DataVersion { $_[0]->_FeatClass->DataVersion }


### "protected"

our $_FeatClass = 'Bio::Chromo::CDS';

sub _FeatClass { $_FeatClass }

sub _info {
    warn("$_[1]\n")
}

sub _warn {
    require Carp;
    Carp::cluck("!! $_[1]\n")
}

sub _croak {
    require Carp;
    Carp::confess($_[1]);
}

sub _mkdir {
    my ($self, $dir) = @_;
    unless ( -d $dir ) {
        require File::Path;
        eval { File::Path::make_path($dir, { verbose => 1 }) };
    }
    -d $dir
}

=method _fetch_gb()

    $seq = $self->_fetch_gb($id[, $gb[, $dir]])

Fetch sequence with id $id from Genbank. If $gb is also given, this is the 
full path of the file to which we write the result, for reuse. $dir is the 
directory portion of $gb, which may be provided for efficiency.

Returns the resulting L<Bio::Seq> object, or an empty list if fetching 
failed.

=cut

sub _fetch_gb {
    my ($self, $id, $gb, $dir) = @_;
    require Bio::DB::GenBank;
    my $db = new Bio::DB::GenBank;
    my $seq = $db->get_Seq_by_id($id);
    unless ( $seq ) {
        $self->_warn("Failed to retrieve sequence!");
        return ();
    }
    return $seq unless $gb;
    # write gb file
    (undef, $dir, undef) = splitpath($gb) unless $dir;
    if ( $self->_mkdir($dir) ) {
        $self->_info("  Writing genbank file '$gb'");
        require Bio::SeqIO;
        my $out = new Bio::SeqIO -format => 'genbank', -file => ">$gb";
        if ( $out ) {
            $self->_warn("  Failed to write sequence to '$gb'")
                unless $out->write_seq($seq);
        } else {
            $self->_warn("  Failed to open '$gb' for writing");
        }
    } else {
        $self->_warn(
          "  Failed to create directory '$dir', no genbank files will be stored"
      );
    }
    $seq
}

1; # End of Bio::Chromo
