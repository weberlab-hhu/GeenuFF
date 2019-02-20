# GeenuFF
Schema and API for a relational db that encodes gene models in an explicit, structured, and robust fashion.

## beta disclaimer

GeenuFF is currently _extremely_ beta and very unstable. 
We're keen to get feed back or ideas from the community
(even if it's just whether you think this could be useful to you
if developed further), but if you build on GeenuFF as it is now, 
you're doing so at your own risk.

## What


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
protein (not 'gene') has the protein ID, the gene, which can 
derive from multiple loci (real or artificial) doesn't have
direct assignment of coordinates. 

#### Consistent
naming, encoding / positioning, 

same code parses EUK/PROK/TRANS-SPLICE

##### Coordinate system

#### Extensible
relational db allows for structured extensions (e.g. to add
info for train/dev/test sets when using for machine learning).

## Install
GeenuFF has been tested exclusively (so far) in python3.6

I would recommend installation in a virtual envinronment.
https://docs.python-guide.org/dev/virtualenvs/

The dustdas dependency is available from github.
From a directory of your choice (and preferably in a virtualenv):

```bash
git clone https://github.com/janinamass/dustdas.git
cd dustdas
python setup.py install
cd ..
```

Afterwards install GeenuFF in the same fashion.

```bash
git clone https://github.com/weberlab-hhu/GeenuFF.git
cd GeenuFF
pip install -r requirements.txt
python setup.py install
cd ..
```

And you might want to run the tests (sorry for the strict directory, will fix)
```bash
cd GeenuFF/geenuff
py.test
cd ../..
```
 

