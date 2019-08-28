import os
import pytest
from sqlalchemy import create_engine
from sqlalchemy import func
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import IntegrityError
import sqlalchemy

from .. import orm
from .. import types
from .. import helpers
from ..base.orm import (Genome, Feature, Coordinate, Transcript, TranscriptPiece, SuperLocus,
                        Protein)
from ..base.handlers import SuperLocusHandlerBase, TranscriptHandlerBase
from ..applications.importer import ImportController, InsertCounterHolder, OrganizedGFFEntries
from ..applications.exporters.sequence import FastaExportController
from ..applications.exporters.lengths import LengthExportController
from ..applications.exporter import MODES
from .test_geenuff import mk_memory_session
from ..applications.exporter import TranscriptHandlerBase

EXPORTING_DB = 'testdata/exporting.sqlite3'


@pytest.fixture(scope="module", autouse=True)
def prepare_and_cleanup():
    if not os.getcwd().endswith('GeenuFF/geenuff'):
        pytest.exit('Tests need to be run from GeenuFF/geenuff directory')

    if not os.path.exists(EXPORTING_DB):
        controller = ImportController(database_path='sqlite:///' + EXPORTING_DB)
        controller.add_genome('testdata/exporting.fa', 'testdata/exporting.gff3', clean_gff=True)
    yield
    os.remove(EXPORTING_DB)


def seq_len_controllers():
    econtroller = FastaExportController(db_path_in='sqlite:///' + EXPORTING_DB)
    econtroller.prep_ranges(range_function=MODES['introns'], genomes=None, exclude=None)

    lcontroller = LengthExportController(db_path_in='sqlite:///' + EXPORTING_DB)
    lcontroller.prep_ranges(range_function=MODES['introns'], genomes=None, exclude=None)
    return econtroller, lcontroller


def test_get_intron_seqs():
    """checks that expected intron sequences are produced"""
    # test data has been simplified from an augustus run that previously resulted in erroneous masks
    # and from a partial gene model in the Rcommumnis genome
    econtroller, lcontroller = seq_len_controllers()
    # expect [(fa_id (samtools faidx), length, sequence)]
    # and the fa_id is just for knowing where it went wrong if it does
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 11
    expect = [("Chr1:195000-199000:1543-2114", 572,
               "GTACTTCCAAATCTTCAATTTTGATTCTAAAGATTGGTCCTTTTACTCTGTTTCTCAATT"
               "TGAGTTTTAGGTATTCTTTGATTTTGTATTGGTTTCATTCTAAATATTCATCCTTTACTC"
               "AACTTCTAGATAAGGGATTTAGGTATTCTCAAATTTCCGATTTGATTCCTTTACTCGTTT"
               "CTAGATTGGGGTTTTAGGAATTACCAGTTGGGGGTTTTGCAATTTGCGTAATCAAAGAAT"
               "TTTATTTGTTGTATTGCTTGGTATTGAAGTTTGTCTCTGTTTCTCTACCTCGTCATGTAA"
               "TGTGCTTAGATCCATTAAGTAAATGCTTGTGGATATTTATGTAGATGGTTAAGAGTGATC"
               "GTGATCAGAGTCCTTCTCTTATTTAACTGCATTGCCTGTGAGTTGTGGTCCTGAAGGTTG"
               "TTGTTATTATTGAATTCTATGTATGTATAGATTATGTCATTGGTCTCATGTGGTTTTTAT"
               "GGGTAACGTCTTTACTAATAATAGCACTATGCTTCTGGATTTTGATCTATGTGATCTGTA"
               "ACATTTCTAGTTGGTGTGTCTTTGATTGCCAG"),
              ("Chr1:195000-199000:2220-2299", 80,
               "GTATATATACCGCTGCTCGTATCTCTTTTCCGGTGTTACAAAAGCGATGTCGTGACCTAA"
               "TGCTGGGTTCGTTACTATAG"),
              ("Chr1:195000-199000:2423-2500", 78,
               "GTAAGTCTGGAATAGCTTTTGAGTTGTCCTCTATGTTTATAAGCTATTGTTGTGTGTAAA"
               "CCTTTGTTATATCTGTAG"),
              ("Chr1:195000-199000:2672-2774", 103,
               "GTAAACTATTAAACTCATTAACTCTCTCCTGCAATCTGCAAGGCAGTCTTTAGGAATGTG"
               "AATATTAGGAAATAACTTTTACTTTGTGGGTTGATTTGTTTAG"),
              ("Chr1:195000-199000:2905-2973", 69,
               "GTATTTGATAAGCATATTCAATCCTCTCTATGTTTTTGTAAATGGTCTTACTATGTGTCG"
               "TGGCTGTAG"),
              ("Chr1:195000-199000:3184-3266", 83,
               "GTATACAAAACAATTTGCCTTTACGTTTTTACATTTCTTTAAGAGTTTGAAACATGTCTA"
               "AAGCTGGGATAATATTTTTGCAG")]
    expect += expect[:4] + expect[5:6]
    print(len(expect), 'len expect')
    iexpect = iter(expect)
    ilengths = iter(lcontroller.export_ranges)
    for i, grp in enumerate(econtroller.export_ranges):
        exp = next(iexpect)
        lgrp = next(ilengths)
        seq = econtroller.get_seq(grp)
        length = lcontroller.get_length(lgrp)
        try:
            assert ''.join(seq) == exp[2]
            assert length == exp[1]
        except AssertionError as e:
            print(exp[0])
            print(i)
            raise e


def test_get_exons():
    pass

