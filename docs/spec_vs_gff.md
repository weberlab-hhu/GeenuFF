## Specific changes vs gff
### features

#### types
Feature types differ somewhat from a gff, "geenuff_" has been 
appended to names to reduce confusion where the meaning is not
identical.

##### core types:
* geenuff_transcript: range of pre-mRNA / unspliced transcript
* geenuff_cds: range between the start and stop codon (ignoring introns)
* geenuff_intron / trans intron: range between a donor and acceptor splice site

##### error types:
These do not exist in a gff, but are used in geenuff to denote things
that might be ambiguous or unknown about a gene model. 

Currently errors are assigned when any obvious gene model inconsistency
is encountered during gff parsing. Most error types
are extended from the end of a known feature half way to the next
gene model, while some internal errors (e.g. too_short_intron) can
be assigned more precisesly. If gff format were not used as an intermediary
and gene annotation was performed and stored directly in a geenuff structured
database, all the errors could be assigned more precise ranges for any
ambiguity.

Types are:

* missing_utr_5p
* missing_utr_3p
* empty_super_locus
* missing_start_codon
* missing_stop_codon
* wrong_starting_phase 
* mismatched_ending_phase
* overlapping_exons
* too_short_intron
* super_loci_overlap_error

##### start_is_biological_start and end_is_biological_end:
When `True`, these attributes mean the start and end attributes
of a feature correspond to a meaningful biological transition.
* start: start of a region, inclusive
  * geenuff_transcript --> transcription start site (1st transcribed bp)
  * geenuff_cds --> start codon, the A of the ATG
  * geenuf_intron --> donor splice site, first bp of intron
* end: end of a region, exclusive (i.e. start of one there after)
  * geenuff_transcript --> 1 after transcription termination site (1st non-transcribed bp)
  * geenuf_cds --> 1 after stop codon, first non-coding bp, e.g. the N in TGAN
  * geenuff_intron --> 1 after acceptor splice site (1st bp that is part of final transcript)

When `False`, these attributes mean the start and end attributes
of a feature either do not, or it is not known if they correspond
to a biological transition, yet the region they delineate is
confidently of the given type. 

For instance, if the parser finds a gene model in a gff where the
start of the first exon and the start of the first CDS (+ strand)
have the same position (the A in ATG), then it is apparent that
we are missing the 5' UTR, so for the geenuff_transcript feature
the start_is_biological_start will be set to False, and an
error mask will be added upstream of the CDS. We are still confident
that all of the CDS must occur within the transcript, we know
the start codon is part of the transcript region, but we mark that the start point
itself is probably wrong, and mask the upstream range as it's 
unclear what part of this is intergenic and which part UTR.
 
#### feature start/end/at numbering

Features have start and end coordinates that
to delineate a range. 

The positioning of these features is in keeping with the common
coordinate system: count from 0, start inclusive, end exclusive. 
So, the "geenuff_cds, start", is at the A, of the ATG, AKA the first
coding base pair; while in contrast, the "geenuff_cds, end" is
after the stop-codon, AKA, the first non-coding bp.

##### reverse complement

Importantly, the coding-start should always point to the first
A, of ATG, regardless of strand. This means the numeric coordinates
have to change and unfortunately while one could take the
sequence \[1, 4) on the + strand, and directly use 1 and 4 as python coordinates
and get the sequence; the same is not going to work on the minus strand.
Instead: 

```
 0  1  2  3  4  5
.N [A .T .G )N .N
 |  |  |  |  |  |
 N. T. A. C. N. N.
```

To get the reverse complement of this on the minus strand, we set the
inclusive start to 3, and exclusive end to 0. Note this is now off by 
one from the python coordinates


```
 0  1  2  3  4  5
.N [A .T .G )N .N
 |  |  |  |  |  |
 N( T. A. C] N. N.
```

##### differences vs gff
Cheat sheet for how the Features compare to the gff (in particular any discrepancy
between the closest coordinate in the gff, and the now standardized, consistent coordinate).

First and last for gff are reported as they are typically in gff (coordinate sorted),
so reverse to the interpretation when on the - strand. 

Plus strand (+)

| Common Name  | GFF | GFF start | GFF end |geenuff type| bearing| position |
| -------------|:----| ---------:|--------:|:---|:-------|--------:|
| TSS, Transcription start site      | start 1st exon  |x| |transcript |start|x - 1|
| TTS, Transcription termination site| end last exon   | |x|transcript |end  |x    |
| 1st bp of start codon              | start 1st CDS   |x| |cds        |start|x - 1|
| coding end                         | end last CDS    | |x|cds        |end  |x    |
| donor splice site (5' of intron)   |end non-last exon| |x|intron     |start|x    |
| acceptor splice site (3' of intron)|start 2nd+ exon  |x| |intron     |end  |x - 1|


Minus strand (-)

| Common Name  | GFF | GFF start | GFF end |genuff type| bearing| position|
| -------------|:----| ---------:|--------:|:---|:-------|--------:|
| TSS, Transcription start site      | end last exon   | |x|transcript |start|x - 1|
| TTS, Transcription termination site| start 1st exon  |x| |transcript |end  |x - 2|
| 1st bp of start codon              | end last CDS    | |x|cds        |start|x - 1|
| coding end                         | start 1st CDS   |x| |cds        |end  |x - 2|
| donor splice site (5' of intron)   |start 2nd+ exon  |x| |intron     |start|x - 2|
| acceptor splice site (3' of intron)|end non-last exon| |x|intron     |end  |x - 1|

