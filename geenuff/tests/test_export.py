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

from .test_geenuff import mk_memory_session
from ..applications.exporter import TranscriptHandlerBase


def test_get_intron_seqs():
    """checks that expected intron sequences are produced"""
    # test data has been simplified from an augustus run that previously resulted in erroneous masks
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/exporting.fa', 'testdata/exporting.gff3', clean_gff=True)