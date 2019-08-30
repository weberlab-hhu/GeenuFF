# GeenuFF

## Source & Install
See [github repository](https://github.com/weberlab-hhu/GeenuFF)

## What
Relational db Schema & Api to store and interpret gene structure

### Conceptual summary
GeenuFF format essentially defines ranges on the genomic sequence 
of a gene structure related "type" (e.g. transcribed, coding, ...).
These ranges are delineated by a feature with a type, start and end. 
The combination of features ultimately required to define
a processed biological macromolecule (e.g. mRNA, protein) from the genomic
sequence, is recorded in links to the 'outer' tables (e.g. "transcript", 
"protein"). 

GeenuFF is designed to unambiguously encode even 
partial information. For instance, the "start_is_biological_start"
and "end_is_biological_end"
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
via trans-splicing that truly derives from two (or more) loci.

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
adding one more standard to the [pile](https://xkcd.com/927/),
then the following section explains _why_ we think it's worth it.

Historically, gene annotations have been encoded in some variation
of the gff format (gtf, gff, gff3), see
[http://gmod.org/wiki/GFF3](http://gmod.org/wiki/GFF3).
There are a variety of drawbacks to these formats. Some are somewhat
more superficial, such as the custom encoding of key, value
pairs specifying both relationships and extra meta info 
in the final column. A good base parser, or non-text alternatives 
such as gffutils
([https://daler.github.io/gffutils/)](https://daler.github.io/gffutils/));
are sufficient
to address such issues. Here, we will focus on the benefits of
the basic changes in underlying _structure_ that address the more
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
in one of the _implicit_ features, the greater problem occurs
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
While GeenuFF also ecodes ranges, and to avoid redundancy still
has some _implicit_ encodings, the choice of _which_ structural
elements to encode has been done to better reflect biology
and so that the start and end of features have a consistent interpretation
for each type. 

For instance, the same two-exon, + strand transcript
used above would basically change as follows
```
# gff-like
 [      transcript        ]
 [ exon ]          [ exon ]

# geenuff 
 [       transcript        )
         [ intron  )
 ^
 start, TSS
```
Or on the - strand:
```
# geenuff 
(       transcript        ]
        ( intron  ]
                          ^
                          start, TSS
```
This already makes it easier to parse, e.g. you don't 
have to use a different rule to find the transcription start site
on the minus strand. 

More importantly however, both start and end, come with an additional
boolean attribute (start_is_biological_start and end_is_biological_end)
which indicate whether this transition is a biological one (when
set to `True` or whether we have a partially known gene model,
when set to `False`).

If we now return to how to encode our incomplete gene model, where
we know we are missing the first exon, we can do so as follows.
```
# geenuff
                       [     transcript      )
                       ^                     ^
                       start                 end

start_is_biological_start=False
end_is_biological_end=True

# further, if we want to indicate our lack of knowledge on the region
# before (perhaps we don't know for sure if this is an intron or an assembly
# error); we can add an error mask feature (by choosing our specific error as "type") 
to indicate our uncertainty. e.g. 
[    missing_utr_5p    )
^                      ^
start                  end
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
field (Derives_from=<a gene ID>) to determine which gene they are associated
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
* it's very unclear whether a protein derived from trans-splicing is part 
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
  (directly or indirectly) to all of the following:
  * Feature: this holds things like "geenuff_transcript", "geenuff_cds", or "geenuff_intron"
    that can directly be assigned coordinates on the genome.
    * appending "geenuff_" is done to disambiguate them from the similarly named gff features
  * Transcript: (AKA pre-mRNA), has an ID, importantly it links 
    (via intermediary, see [transcript_piece](#transcript_piece))) 
    to all the features combined to make a transcript and 
    any protein translated there from.
  * Protein: (AKA protein), has an ID, importantly it links to geenuff_cds-type features 
    making one protein.

The gains of this restructuring is that the gene structure of Eukaryotic, Prokaryotic,
and even Trans-spliced examples can be encoded in fundamentally the same fashion and
parsed with the same code.

For example:
* In a common Eukaryotic example, a SuperLocus will point to multiple Transccripts. 
Each Transcript will have one "geenuff_transcript" and one "geenuff_cds" feature
and will be connected to one Protein.

* In a common Prokaryotic example, a SuperLocus will point to one Transcript
with one "geenuff_transcript" feature, but
this Transcript will be connected to multiple Proteins and have
multiple "geenuff_cds" features.

* In a case of trans-splicing, a Transcript will have "geenuff_transcript" features,
one at it's different loci. How these 
pieces are ultimately linked is encoded using the "position" attribute of the
TranscriptPieces (see [transcript_piece](#transcript_piece)
and [feature](#feature)).

* Finally, any of the above could be encoded across multiple artificial breaks in the
sequence (such as a fragmented assembly), by setting the "<>\_is_biological\_<>" 
attributes, adding error masks as necessary and using TranscriptPieces 
(see [transcribed_piece](#transcribed_piece) and [feature](#feature))).

While these four examples differ in which features they use, they all follow the
same spec and can be parsed with the same code. For instance, ("geenuff_cds", start) is always
interpreted in the exact same way; and you don't get cases like those in a gff, where
the edge of a CDS feature means something else depending upon the presence of other
'CDS' features and the relative position of the overlapping 'exon' feature.
Similarly, the same code can be used to interpret a split-locus gene model, whether
this has a biological (trans-splicing) origin or artificial (e.g.
fragmented assembly) origin. The features differ as necessary, the logic remains the same.

#### gff-like coordinate troubles.
The most common usage of a gff-like file is simply to denote which sequences
belong to the original transcript, the final transcript and the proteins. This
works decently with the `[inclusive start, inclusive end]` coordinates the gff-like
formats use. However, the _implicit_ components are then inverted, and tedious, e.g.
for an intron the coordinates are `(exclude end exon_i, exclude start exon_i+1)`,
or for the 3' untranslated region (UTR) the coordinates are 
`[inclusive start exon, exclude start CDS)`.

##### GeenuFF increased coordinate consistency
In the geenuff format the start of features are always inclusive
while the end of features are always exclusive. Besides being more
consistent with e.g. python coordinates, this has the major advantage that if one
wants a component that isn't explicitly included for reasons of avoiding redundancy
(e.g to get the UTR you still have to take "geenuff_transcript" - "geenuff_cds") 
it's at least always 
`[inclusive start, exclusive end <or start next>)`. So the first coding exon 
would be `[geenuff_cds start, geenuff_intron start)`. 

__Caveat:__ 
Ranges on the minus strand are off by one from the pythonic coordinates.

#### Extensible
Finally, with geenuff being based on a relational database, it's much easier
to _extend_ the format for specific purposes without changing the base, shared
attributes. 

For instance, we needed to track some meta_information
for an applied machine learning and gene model project. This required just
an additional tables with a foreign key to Coordinate, but required no modification
of the core format.

## What (with details)

### Comparison of spec to gff
Quick start / reference for 1:1 comparison with gff:
[spec_vs_gff.html](spec_vs_gff.html)

### Coordinate system
GeenuFF coordinates count from 0, have an inclusive start and exclusive end.

The position of the _start_ attribute always marks the 
_inclusive_ start of this range. The position of the _end_ attribute
always marks the _exclusive_ end of the range. 

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
genomic_sequence[end + 1, start + 1]
```
of course one will want to reverse and complement
this sequence as well for most purposes, so a handling function is anyways advisable. 
(and available in `applications.exporters.sequence`)

### Schema summary
For the fine details, please look at `base/orm.py`; here we will just
try and describe the overall structure / major pieces, and what they mean.

#### Tables & relations
Indentation inside a piece indicates a one-to-? relation

* genome
  * coordinate
* super_locus
  * transcript, < many2many to protein >
    * transcript_piece, < many2many to feature >
  * protein, < many2many with transcript, feature >
  * feature, < many2many with protein, transcript_piece; many2one to coordinate >

(and linkage-only association tables for the many2many fields)

##### genome
This _mostly_ holds meta information 

##### coordinate
(sequence meta info: seqid, length
sha1 hash of the sequence the annotation is for,
optional full sequence)

#### super_loci and children
###### super_loci
essentially delineates the graph of things that might possibly be combined,

has given_name and type for the ~gene

###### feature
these describe geenuff specific types and ranges, but
otherwise resemble features from a gff

* they occur on a sequence (foreign key to coordinates)
* they have a position (start and end)
* they have a type (or channel), e.g. {geenuff_transcript, geenuff_cds, geenuff_intron, various error types...}
* they have a start_is_biological_start and end_is_biological_end which indicates whether 
start and end mark the biologically meaningful transition or just the edge of what we know.
* they have a boolean indicator is_plus_strand.
* they have a score (confidence)
* they have a phase (important for type "geenuff_cds", else None; usage and checking still need to be implemented, todo). This is particularly
important if trying to encode a partial protein sequence.
* they have a given_name (but this is often "None" as nothing was available)
* they have a source
* they have many to many relationship with protein (mostly to assign the 'protein_id')
* they have a many to many relationship with transcript_piece.
at the first bp of the downstream sequence.

__Feature self consistency:__
All features on a transcribed piece must be interpretable when sorted 5' to 3'.
There are some biological considerations
such as: an intronic (_cis_ or _trans_) range cannot have a "start" or "end" 
unless it's in a transcribed region, a 
coding range can only "start" or "end" inside a transcript, but non-intronic region.
All features assigned to a transcribed_piece must have the same value for "is_plus_strand".
If a feature has a False value for <>\_is_biological\_<>  attribute, these should occur at the edge
of the piece, or be accompanied by a feature of with an error type to mask the ambiguous area.

###### protein
Each Protein object basically just points to one protein's worth of "geenuff_cds" 
type features (and associated transcript & super_locus),
and has a given_name attribute.

###### transcript_piece
TranscriptPieces delineate a collection of features that can be interpreted 
(5'-3') together. In the standard case (with either _cis_- or no- splicing, 
and representing a full gene model)
a Transcript will have a single TranscribedPiece pointing to all the 
relevant Features. 

In cases where the Transcribpt (final mRNA) is split for either biological or 
technical reasons, each involved locus should have its own transcript_piece.

Each TranscriptPiece should have a single "geenuff_transcript"
feature covering its whole range;
that is the most 5' part of the piece should be at the _start_ from the
"geenuff_transcript" feature and the most 3' part of the piece at the _end_.
The exception is error type features, which may be associated with the piece,
but excede this range.

transcribed_piece has a many2one relationship with transcript.

The 5' to 3' ordering of transcript_pieces within a transcript can be accomplished using
the 'position' attribute of the transcript_piece.

###### transcript
Transcript objects ultimately define what should be interpreted together
to produce the final biological molecule (e.g. pre-mRNA, mRNA, protein, etc..).

They consist of one or more transcript_pieces (which can be ordered 5' to 3' by sorting 'position'
in ascending order). The features
within each transcribed_piece can simply be ordered by coordinates. With pieces
and their features sorted, a transcript can be read 5'-3' and the information
of interest (beit the whole transcript range, the spice sites, the start codon, etc..)
can be extracted as necessary. 

Example logic for interpreting a transcript is can be found in
in `geenuff.applications.exporter.RangeMaker`.

