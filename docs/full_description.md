# GeenuFF

## What
Relational db Schema & Api to store and interpret gene structure

### Conceptual summary
GeenuFF format essentially defines ranges on the genomic sequence 
of a gene structure related "type" (e.g. transcribed, coding, ...).
These ranges are delineated by one each of a starting & ending feature
of their type. The combination of features ultimately required to define
a processed biological macromolecule (e.g. mRNA, protein) from the genomic
sequence, is recorded in the 'outer' tables (e.g. transcribeds, 
translateds) that link to features. 

GeenuFF is designed to unambiguously encode even 
partial information. For instance, the "bearing"
of a feature can be used to indicate whether the feature delineates
the expected full biological range, or merely a part of it. Thus a 
gene model encoded in geenuff can explicitly differentiate a full
from partial gene model, and exactly _where_ our knowledge of the
gene model ends.

GeenuFF is designed to encode all necessary complexity, both to 
represent the biological processes, as well as technical limitations.
In particular, the abstraction layer of transcribed_pieces can be used
to group and organize features
encoding a transcript that originates from multiple unrelated genomic loci.
Thus geenuff can encode a transcript split across two scaffolds in a 
fragmented genome assembly, or it can encode a transcript processed
via trans-splicing, which has a true biological origin from two (or more) loci.

### Coordinate system
GeenuFF coordinates count from 0, have an inclusive start and exclusive end.

GeenuFF features always come in pairs defining a range. The position of the
_starting_ feature (bearing in ["start", "status_open"]) always marks the 
_inclusive_ start of this range. The position of the _ending_ feature
(bearing in ["end", "status_close"]) always marks the _exclusive_ end
of the range. 

Therefore, GeenuFF coordinates exactly match (among other languages) 
python coordinates for ranges on the positive strand.

However, the logic of inclusive start, exclusive end applies even when encoding
features on the _minus_ strand.

Let's say you want to encode the range to select the elements `{2, 3}` 
(of `[0, 1, 2, 3, 4, ...]`) with geenuff features / coordinates. 
The plus strand would exactly match the pythonic coordinates `[2, 4)`,
but the negative strand (want `{3, 2}`) would include the start point,
and not the end `[3, 1)`. Thus, if the genomic sequence was represented
in a python list you would select a range on the minus strand 
with something like this:
```
genomic_sequence[end_feature.position + 1, start_feature.position + 1]
```
of course one will want to reverse and complement
this sequence as well for most purposes, so a handling function is anyways advisable. 

### Schema summary
For the fine details, please look at `base/orm.py`; here we will just
try and describe the overall structure / major pieces, and what they mean.

#### Tables & relations
Indentation inside a piece indicates a one-to-? relation

* annotated_genomes
  * sequence_infos
    * coordinates
* super_loci
  * super_loci_aliases
  * transcribeds, < many2many to translateds >
    * transcribed_pieces, < many2many to features >
  * translateds, < many2many with transcribeds, features >
  * features < many2many with translateds, transcribed_pieces; many2one to coordinates >
    * upstream_features < via up_down_pairs to downstream_features >
    * downstream_features < via up_down_pairs to upstream_features >
* up_down_pairs < links upstream_features, downstream_features >

(and 3 more linkage-only association tables for the many2many fields)

#### annotated_genomes & children
These _mostly_ hold meta information 
###### annotated_genomes 
(meta info only: species, accesion, version, acquired_from)
###### sequence_infos 
(an abstraction layer, for instance if you wanted to split
within genome sequences for any processing reasons; e.g. denote what's nuclear
vs organellar...)
###### coordinates 
(sequence meta info: seqid, inclusive start, exclusive end)

_Note: should probably add a hash of the sequence or something here too._

#### super_loci and children
###### super_loci
essentially delineates the graph of things that might possibly be combined,

has given_id as meta_info for the ~gene

###### features
these describe transitions that can occur in a gene model but
otherwise resemble features from a gff

* they occur on a sequence (foreign key to coordinates)
* they have a position (currently called start, todo...)
* they have a type (or channel), one of {transcribed, coding, intron, trans_intron, error}
* they have a bearing {start, status_open, status_close, end} which indicates whether something
begins {start, status_open} or ends here {end, status_close}, and whether it is the biologically
meaningful transition {start, end} or should just be used to set the status / is handling partial
information {status_open, status_close}.
* they have a boolean indicator is_plus_strand.
* they have a score (confidence)
* they have a phase (calculation and usage still need to be implemented, todo). This is particularly
important if trying to encode a partial protein sequence.
* they have a given_id 
* they have a source
* they have two additional subtypes: upstream and downstream features (more info below)
* they have many to many relationship with translateds (mostly to assign the 'protein_id')
* they have a many to may relationship with transcribed_pieces.

__Feature pairing:__ 
Features _must_ always delineate a range / come in pairs on a 
trabscribed_piece. So, where there is a feature of type X with bearing "status_open" or 
"start", then there must be a feature downstream (3') on the same transcribed_piece,
of this that is also of type X and is of either bearing "status_close" or "end". 
Similarly, there may be no feature with the bearing "end" or "status_close" without
a corresponding feature with the bearing "status_open" or "start" upstream (5') of
it.  

__Upstream & Downstream Features:__
These feature subtypes are used to order transcribed_pieces within a transcript. If
a transcript has multiple transcribed_pieces (e.g. in the case of trans-splicing), 
two things need to be tended to. First, by the pairing rule above, any status type
that _started_ on the first piece must have a paired _ending_ feature on the same 
piece; while on the second piece, any _ending_ feature must have the accompanying 
_starting_ feature. This means that at the end of the first piece (e.g. the transcription 
termination site) all other 'open types' (at least "trans_intron" for trans-splicing)
must be closed, (of bearing "status_close". This closing feature gets the subtype 
UpstreamFeature, and is linked via the UpDownPair table to a DownstreamFeature. 
The linked DownstreamFeature has the bearing "status_open" and is located at the 
start of the next piece (e.g. transcription start site). While a given transcribed_piece
can have multiple UpstreamFeatures or DownstreamFeatures, it cannot have conflicting 
information (i.e. all upstream features must be at the same position, and point to a
downstream features that are on the same position of the same piece). Finally,
circular links aren't allowed (if Piece A is upstream of Piece B, then Piece B may not,
directly or indirectly, be upstream of Piece A).

Below is an example of how a non-coding, all + strand, trans-spliced transcript would be encoded.
```
# where '~' marks parts of the final transcript, and '.' marks what's spliced out
Piece A:                                         Piece B:
5'                          3'                 5'                           3'
~~~~~~~~....................                   ...................~~~~~~~~~~
<transcribed>
^                           ^                  ^                            ^
start                       end                start                        end

<trans_intron>
        ^                   ^                  ^                  ^  
        start               status_close  <->  status_open        end

```

__Feature reading:__
All features on a transcribed piece must be interpretable when sorted 5' to 3'.
That includes the pairing mentioned above, as well as some biological considerations
such as: an intronic (_cis_ or _trans_) range cannot "start" or "end" 
unless it's in a transcribed region, a 
coding range can only "start" or "end" inside a transcribed, but non intronic region.
All features assigned to a transcribed_piece must have the same value for "is_plus_strand".
Status marking features should occur at the edges of the piece, or be accompanied
by a range of type "error" to mask the ambiguous area.

###### translateds
Each Translated object basically just points to one protein's worth of "coding" 
type features (and associated transcribed & super_locus),
and has a given_id attribute.

###### transcribed_pieces
TranscribedPieces delineate a collection of features that can be interpreted 
(5'-3') together. In the standard case (_cis_- or no- splicing, full gene model)
a Transcribed object with have a single TranscribedPiece pointing to all the 
relevant Features. 

In cases where the Transcribed (final mRNA) is split for either biological or 
technical reasons, each involved locus should have its own transcribed_piece.

The entire region of a TranscribedPiece should be interpreted as "transcribed";
that is, the most 5' piece should be of type "transcribed" and bearing
"start" or "status_open", while the most 3' piece should be of 
type "transcribed" and bearing "end"
or "status_closed" (ties are OK and expected for the status bearings).

See above Feature section for further rules on features of a transcribed_piece.

transcribed_pieces have a many2one relationship with transcribeds

###### transcribeds
Transcribed objects ultimately define what should be interpreted together
to produce the final biological molecule (e.g. pre-mRNA, mRNA, protein, etc..).

They consist of one or more transcribed_pieces (which can be ordered according)
to any Upstream / Downstream Features linked by an UpDownPair. The features
within each transcribed_piece can simply be ordered by coordinates. With pieces
and their features sorted, a transcript can be read 5'-3' and the information
of interest (beit the whole transcript range, the spice sites, the start codon, etc..)
can be extracted as necessary. 

The basic logic for interpreting a transcript is coordinated 
in `geenuff.base.api.TranscriptInterpBase`.

## Why
### general goal
The encoding of gene annotations (which parts of an organisms
genome will be transcribed into RNA, how the RNA will be modified
to make the final mature RNA (mRNA), and what part of the mRNA
will be transcribed into protein) ultimately requires a 
non-trivial and tailored data structure to properly capture the 
information and recreate at will the coordinates of a gene, 
gene part, or the final or intermediate sequences created, etc.
Some of the challenges include handling cases like the following
* where there are multiple ways to 'put together' the transcripts or proteins
originating from a single loci
* multiple transcripts can produce the same final protein and differ only in mRNA
* one transcript (in prokaryotes) can be translated into several proteins
* in rare cases of trans-splicing a single mRNA can be derived from multiple
  genomic loci
* and all of the above have to be encodable in cases of partial information, including:
  * incomplete genomic sequence
  * fragmented genomic sequence
  * incomplete information about the gene annotation itself

### specific advantages over alternatives
Unsurprisingly, GeenuFF does not represent the first attempt
to encode gene structures; so if you're already familiar 
with e.g. the gff format and are wondering if we're just 
adding one more standard to the pile (https://xkcd.com/927/), 
then the following section explains _why_ we think it's worth it.

Historically, gene annotations have been encoded in some variation
of the gff format (gtf, gff, gff3), see http://gmod.org/wiki/GFF3. 
There are a variety of drawbacks to these formats. Some are somewhat
more superficial, such as the custom encoding of key, value
pairs specifying both relationships and extra meta info 
in the final column. A good base parser, or non-text alternatives 
such as gffutils (https://daler.github.io/gffutils/); are sufficient
to address such issues. Here, we will focus on the benefits of
the basic changes in encoding _structure_ that address the more
fundamental issues.

#### Gff-like implicit encodings
In gff-like formats, genes and gene pieces are encoded as 
ranges, e.g. `[ exon ]`. Unfortunately this leaves many gene 
components to only be encoded implicitly. For instance,
an intron is encoded as the gap between two exons:
```
# gff features
[      transcript        ]
[ exon ]          [ exon ]
# interpretation
[ exon ]( intron )[ exon ]
```
Similarly, the transcription start site (TSS) is not indicated
explicitly, but is rather implicitly assumed to be:
```
# at the start of the 1st exon (+ strand)
[    transcript    ]
[ exon ]    [ exon ]
^
TSS

# or, at the end of the last exon (- strand)
[    transcript    ]
[ exon ]    [ exon ]
                   ^
                   TSS
```
While this always requires some extra parsing if one is interested
in one of the _implicit_ features; the greater problem occurs
in cases where one has less-than perfect information. 

For example, lets assume we want to encode a partial gene model;
we know from homology comparison to other species and the truncated
mapping of RNAseq reads that we are missing at least the first exon.
With gff-like formats one can either include the exons one knows, 
which will erroneously _imply_ the transcription start site is located
where one actually has an acceptor splice site; or one can skip
the gene model entirely, which will erroneous _imply_ that the whole
region is intergenic. There is no 'right' way to document what one 
knows and what one doesn't with these formats.

##### GeenuFF more explicit encodings
With GeenuFF, the transitions are encoded explicitly, and are strictly
paired to denote the ranges. For instance, the same two-exon, + strand transcript
used above would basically change as follows
```
# gff-like
[      transcript        ]
[ exon ]          [ exon ]

# geenuff 
<transcribed>
[                         )
^                         ^
start                     end
<intron>
        [         )
        ^         ^
        start     end   
```
This already makes it easier to parse, e.g. you don't 
have to use a different rule to find the transcription start site
on the minus strand. 

More importantly however, in addition to the bearings 'start' and
'end', there are also 'status_open' and 'status_close' to denote
incomplete information as well as 'point' (not currently used).

If we now return to how to encode our incomplete gene model, where
we know we are missing the first exon, we can do so as follows.
```
# geenuff
<transcribed>
        [                     )
        ^                     ^
        status_open           end
# further, if we want to indicate our lack of knowledge on the region
# before (perhaps we don't know for sure if this is an intron or an assembly
# error); we can use the error channel (type) to indicate our uncertainty
<error>
[        )
^        ^
start    end
```

#### Gff-like miss assignment of attributes and relations
Ultimately a gff-like encoded gene model tries to map biology
onto a miss-fitting model. 

For instance, in a gff-like model a transcript is always the
child of a gene. This works fine for most Eukaryotic genes, but in prokaryotes
where several proteins are derived from a single transcript; this 
requires an awkward patch to encode in a gff file. Specifically,
the genes / proteins derived from one transcript are labeled
as children of a new feature, the operon. The single transcript 
then has multiple parent genes, and the CDS pieces now have 
to use a prokaryote-specific key-value pair in their attribute
field (Derived_from=<a gene ID>) to determine which gene they are associated
with. Note that this drastically changes the meaning of the _implicit_
features discussed above, and requires fundamentally different code to
parse / interpret.

Even more problematic are 'corner cases' like trans-splicing,
where one final mature mRNA, may be derived from two distant
original transcripts which are then ligated together. From our experience,
we have not seen a standard way, and there certainly isn't a good
way to encode this in a gff-like format.

The miss-fitting relationship structure is further exacerbated by 
assigning attributes at the wrong level. For instance,
trans-splicing makes it clear that it is not the gene nor protein that
should have on-genome coordinates; but that they are made up of pieces,
such as transcription or coding start and end sites, that have on-genome
coordinates. Another good example of the miss-placed assignment is common
(yet not standardized) practice of deriving the protein ID from the gene 
or transcript ID, instead of assigning it specifically to the protein. 

While some of this confusion most certainly derives from a somewhat 
undefined biological concept of a gene...
* 'gene' is frequently used to refer to a genomic locus so that one 'gene ID' can
  be assigned to the often highly-related transcripts 'put together differently'
  by alternative splicing. In such cases, it's an umbrella of sorts that encompasses
  all mRNA _transcribed_ from that locus, as well as all proteins _translated_ from
  the mRNA.
* yet in prokaryotes 'gene' is used more to refer to the proteins and rather denotes
  a sub-section of what is _transcribed_.
* while it's very unclear whether a protein derived from trans-splicing is part 
  of one or two 'genes'?

... be this as it may, this is no excuse to encode what is unambiguous in a
clear and consistent fashion. 


##### GeenuFF restructuring to bring the map closer to the territory
GeenuFF essentially consists of a schema for a relational database (and an API to 
interpret as necessary). The schema breaks up the artificial connection
rules of the gff-like formats; and with many-to-many fields directly allows
assignment of things like multiple transcripts to one protein, or multiple 
proteins to one transcript. The key tables / concepts in the schema are as follows.

* SuperLocus: this is a holder for related transcripts, proteins, and all the bits
  that might be combined to make them. Essentially, it's an artificial / abstract
  concept, but importantly it delineates the maximum graph
  one might have to walk to to put any-subcomponent in context, by linking 
  (one-to-many) to all of the following:
  * Feature: this holds things like "coding, start" (start codon), or "intron, end"
    (acceptor splice site), that can directly be assigned coordinates on the genome.
  * Transcribed: (AKA transcript / mRNA), has an ID, importantly it links 
    (via intermediary, see pieces spec) to all the features combined to make a transcript and 
    any protein translated there from.
  * Translated: (AKA protein), has an ID, importantly it links to coding-type features 
    making one protein.

Note that we have specifically avoided using terms like 'transcript' or 'gene' which
are used in a gff for a related, yet slightly different concept.

The gains of this restructuring is that the gene structure of Eukaryotic, Prokaryotic,
and even Trans-spliced examples can be encoded in fundamentally the same fashion and
parsed with the same code.

For example:
* In a common Eukaryotic example, a SuperLocus will point to multiple Transcribeds. 
Each Transcript will have one protein,
and therefore one each of "coding, start" and "coding, end" Features.

* In a common Prokaryotic example, a SuperLocus will point to one Transcribed, but
this Transcribed with have multiple proteins, and therefore multiple pairs
of "coding, start" and "coding, end" features.

* In a case of trans-splicing, a Transcribed will have multiple paired features for
"transcribed, start" and "transcribed, end" at it's different loci. How these 
pieces are ultimately linked is encoded with the addition of "trans_intron" 
and "status" up/downstream features (see pieces spec).

* Finally, any of the above could be encoded across multiple artificial breaks in the
sequence (such as a fragmented assembly), with the addition of the "status" up/downstream
features (see pieces spec).

While these four examples differ in which features they use, they all follow the
same spec and can be parsed with the same code. For instance, "coding, start" is always
interpreted in the exact same way; and you don't get cases like those in a gff, where
the edge of a CDS feature means something else depending upon the presence of other
'CDS' features and the relative position of the overlapping 'exon' feature.
Similarly, the same code can be used to interpret a split-locus gene model, whether
this has a biological (trans-splicing) origin or artificial (intentional data split,
fragmented assembly) origin. The features differ as necessary, the logic remains the same.

#### gff-like coordinate troubles.
The most common usage of a gff-like file is simply to denote which sequences
belong to the original transcript, the final transcript and the proteins. This
works decently with the `[inclusive start, inclusive end]` coordinates the gff-like
formats use. However, the _implicit_ components are then inverted, and tedious, e.g.
for an intron the coordinates are `(exclude end exon_i, exclude start exon_i+1`,
or for the 3' untranslated region (UTR) the coordinates are 
`[inclusive start exon, exclude start CDS)`.

##### GeenuFF increased coordinate consistency
In the geenuff format the start (& status_open) features are always inclusive
while the end (& status_close) features are always exclusive. Besides being more
consistent with e.g. python coordinates, this has the major advantage that if one
wants a component that isn't explicitly included for reasons of avoiding redundancy
(e.g to get the UTR you still have to take transcribed - translated in geenuff) 
it's at least always 
`[inclusive start, exclusive end <or start next>)`. So the first exon 
would be `[coding start, intron start)`. More generally, while transitioning
5'-3' one can always take the range indicated by `[feature_i.position, feature_i+1.position)`
and assign to this the status (coding/intronic/UTR/etc...) that was established at feature_i.


#### Extensible
Finally, with geenuff being based on a relational database, it's much easier
to _extend_ the format for specific purposes without changing the base, shared
attributes. 

For instance, we needed to assign sequence regions to training / dev / and test 
sets for an applied machine learning and gene model project. This required just
an additional table with a foreign key to the sequences, but required no modification
of the core format.
