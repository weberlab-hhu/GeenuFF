# Command line tools

__Warning! little about this is stable or tested yet__

## import a species into the database

You can import a Eukaryotic species genome + annotation into geenuff db (`geenuff.sqlite3`)
and log the import to `geenuff.import.log` as follows:

Capital letters and `<>` indicate what must be user specified.

```

import2geenuff.py --fasta <PATH_TO_GENOME_FASTA_FILE> --gff3 <PATH_TO_GFF3_FILE> \
    --db-path geenuff.sqlite3 --log-file geenuff.import.log --species <SPECIES_NAME>
```

Or with some of the testdata filled in:
```
geenuff_path=<PATH/TO/GeenuFF/>

import2geenuff.py --fasta $geenuff_path/geenuff/testdata/exporter.fa \
    --gff3 $geenuff_path/geenuff/testdata/exporter.gff3 \
    --db-path geenuff.sqlite3 --log-file geenuff.import.log --species dummy
```


If it is not a problem to pre-arrange your files into the following structure,
the `--basedir` option can be used to simplify/structure input and output.

```
# target structure
<BASEDIR>/input/<YOUR_FILE>.fa
<BASEDIR>/input/<YOUR_FILE>.gff3

# simplified import
import2geenuff.py --basedir <BASEDIR> --species <SPECIES_NAME>

# the output files "geenuff.sqlite3" and "import.log" will be written in
# the directory <BASEDIR>/output/
```

Any custom input parameters will overwrite those from `--basedir` if both are
specified.

Annotations for multiple species can be stored in one database by using the
same value for `--db-path`.

Two genomes cannot be imported with the same "species" name, and this
will throw an integrity error. You can either delete the database
or set the `--replace-db` parameter.

## extract fasta sequences from a database

You can dump transcript, cds, (hopefully soon protein) and a variety
of other sequence breakdowns from the database back to fasta format using the
script `GeenuFF/scripts/dump_to_fasta.py`

For instance:

```
# to obtain the CDS (ignoring phase) you can run
python $geenuff_path/scripts/dump_to_fasta.py --db-path-in GENUFF_DB --mode CDS

# to obtain the final mRNA sequence you can run
python $geenuff_path/scripts/dump_to_fasta.py --db-path-in GENUFF_DB --mode mRNA

# to obtain the unspliced transcript you can run
python $geenuff_path/scripts/dump_to_fasta.py --db-path-in GENUFF_DB --mode pre-mRNA
```

## extract sequence lengths from a database

See above, but use `$geenuff_path/scripts/dump_lengthinfo.py`

Or directly obtain summary statistics by specifying `--stats-only`

# GeenuFF API

Sometimes one wants to have a little more flexibility than
already implemented. In such case it may be more useful to use
the python module.

Warning: nothing about this is stable.

## import
If you're looking to change how and where a genome can be imported
(we do this a lot for e.g. test cases), you can call the python functions more
directly. E.g. for some of the test data from `$geenuff_path/geenuff`

```{python}
from geenuff.applications.importer import ImportController

controller = ImportController(database_path='sqlite:///' + EXPORTING_DB)
        controller.add_genome('testdata/exporting.fa', 'testdata/exporting.gff3', clean_gff=True,
                              genome_args={'species': 'dummy'})
```

You can also use an in-memory database `'sqlite:///:memory:'`, which
can be nice for testing or exploring the data / data structure.

At some point we will expand to include a prokaryotic gff importer
and to make a nice api where a user can make small changes to 
accommodate their own non-conforming gff (because we really can't
claim to support them all). But we aren't there yet.

## sequence output
look at `geenuff.applications.exporter` and `geenuff.applications.exporters.sequence`

For instance you might want to use or extend the class `FastaExportController`
if you wanted to do something with the sequence breakdowns besides writing to a fasta file.

e.g.

```{python}
from geenuff.applications.exporters.sequence import FastaExportController

controller = FastaExportController(PATH_TO_GEENUFF_DB)
controller.prep_ranges(mode='pre-mRNA')
for export_group in controller.export_ranges:
    pre_mrna_seq = controller.get_seq(export_group)
    # your code here
    # check for your motif of interest, count kmers, send to custom output, etc...
```

If you need to go beyond the currently available sequence breakdowns,
you can look at `RangeMaker` from `geenuff.applications.exporter`.export_group

## gff3 output
todo

## json output

It's about the most naive implementation possible at the moment
but geenuff can now query an region and return a somewhat flattened
json output.

e.g. were PATH_TO_GEENUFF_DB was imported from the test files:
"geenuff/testdata/exporting.*" as in the second example for 
the `import2geenuff.py` section above.

```
from geenuff.applications.exporters.json import JsonExportController

controller = JsonExportController(PATH_TO_GEENUFF_DB)
json_out = controller.coordinate_range_to_json(species='dummy',
                                               seqid='Chr1:195000-199000',
                                               start=1, end=3900, 
                                               is_plus_strand=True)
```
This should return all overlapping super loci and _all_ flattened
 children there of. The returned json should have (hopefully something very 
similar to, I'm sure there's mistakes) the following format:

```
[{"coordinate_piece": 
    {"id": int, "seqid": str, "sequence": str, "start": int, "end": int},
 "super_loci":
    [{"id": str, 
      "given_name": str,
      "is_fully_contained": bool,
      "overlaps": bool,
      "transcripts": [{"id": str,
                       "given_name": str,
                       "is_fully_contained": bool,
                       "overlaps": bool,
                       "features": [{"id": int,
                                     "given_name": str,
                                     "seqid": str,
                                     "protein_id", str,
                                     "type": str,
                                     "start": int,
                                     "start_is_biological_start": bool,
                                     "end": int,
                                     "end_is_biological_end": bool,
                                     "score": float,
                                     "source": str,
                                     "phase": int (in {0, 1, 2}),
                                     "is_plus_strand": bool,
                                     "is_fully_contained": bool,
                                     "overlaps": bool
                                    }, ...]
                      }, ...]
     }, ...]

}, ...]
```

All structure (e.g. many to many relationships) in the database that cannot 
be captured in the above format will be handled by including the lower in the
hierarchy elements redundantly. So if you query a super locus which is split
across two scaffolds, two of the above "coordinate_pieces" will come back
and each will have the full set of super_locus, transcript and feature elements
except that "is_fully_contained" and "overlaps" will be updated for the particular
coordinate_piece. Similarly, any transcripts sharing features will simply have
the features repeated redundantly for each transcript. 

Features in a transcript will always be reported in 5'-3' order.
