# GeenuFF API

## json output
The api will be able to run a query based on e.g. a sequence range
or superlocus id and return a JSON of (hopefully something very 
similar to) the following format:

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
                                     "type": enum,
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

All structure (e.g. many to many relationthips) in the database that cannot 
be captured in the above format will be handled by including the lower in the
hierarchy elements redundantly. So if you query a super locus which is split
across two scaffolds, two of the above "coordinate_pieces" will come back
and each will have the full set of super_locus, transcript and feature elements
except that "is_fully_contained" and "overlaps" will be updated for the particular
coordinate_piece. Similarly, any transcripts sharing features will simply have
the features repeated redundantly for each transcript. 

Features in a transcript will always be reported in 5'-3' order.
