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
 

