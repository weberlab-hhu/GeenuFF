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


@pytest.fixture(scope="session", autouse=True)
def prepare(request):
    if not os.getcwd().endswith('GeenuFF/geenuff'):
        pytest.exit('Tests need to be run from GeenuFF/geenuff directory')


### Helper functions ###
def mk_memory_session(db_path='sqlite:///:memory:'):
    engine = create_engine(db_path, echo=False)
    orm.Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    return Session()


def setup_data_handler(handler_type, data_type, **kwargs):
    data = data_type(**kwargs)
    handler = handler_type()
    handler.add_data(data)
    return data, handler


def cleaned_commited_features(sess):
    all_features = sess.query(Feature).all()
    allowed_types = [
        types.TRANSCRIPT_FEATURE, types.CODING, types.INTRON, types.TRANS_INTRON, types.ERROR
    ]
    clean_datas = [x for x in all_features if x.type.value in allowed_types]
    return clean_datas


# Do not implement as __eq__ to not change __hash__ behavior
def matching_super_loci(s1, s2):
    if not isinstance(s1, s2.__class__):
        return False
    if s1.id is not None and s2.id is not None:
        if s1.id != s2.id:
            return False
    if s1.given_name != s2.given_name or s1.type != s2.type:
        return False
    return True


def matching_features(f1, f2):
    if not isinstance(f1, f2.__class__):
        return False
    if f1.id is not None and f2.id is not None:
        if f1.id != f2.id:
            return False
    if (f1.type != f2.type or f1.given_name != f2.given_name or f1.start != f2.start
            or f1.end != f2.end or f1.start_is_biological_start != f2.start_is_biological_start
            or f1.end_is_biological_end != f2.end_is_biological_end
            or f1.is_plus_strand != f2.is_plus_strand or f1.phase != f2.phase
            or f1.coordinate.id != f2.coordinate.id):
        return False
    return True


def matching_proteins(t1, t2):
    if not isinstance(t1, t2.__class__):
        return False
    if t1.id is not None and t2.id is not None:
        if t1.id != t2.id:
            return False
    if t1.given_name != t2.given_name or t1.super_locus.id != t2.super_locus.id:
        return False
    return True


def matching_transcripts(t1, t2):
    if not isinstance(t1, t2.__class__):
        return False
    if t1.id is not None and t2.id is not None:
        if t1.id != t2.id:
            return False
    if (t1.given_name != t2.given_name or t1.super_locus.id != t2.super_locus.id
            or t1.type != t2.type):
        return False
    return True


def orm_object_in_list(obj, obj_list):
    for o in obj_list[:]:  # make a copy at each iteration so we avoid weird errors
        if ((isinstance(o, Feature) and matching_features(o, obj))
                or (isinstance(o, Protein) and matching_proteins(o, obj))
                or (isinstance(o, Transcript) and matching_transcripts(o, obj))):
            obj_list.remove(o)
            return True
    return False


### The actual tests ###
def test_annogenome2coordinate_relation():
    """Check if everything is consistent when we add an Genome and a Coordinate
    to the db. Also check for correct deletion behavior.
    """
    sess = mk_memory_session()
    g = Genome(species='Athaliana', version='1.2', acquired_from='Phytozome12')
    coord = Coordinate(seqid='abc', genome=g)
    assert g is coord.genome
    # actually put everything in db
    sess.add(coord)
    sess.commit()
    # check primary keys were assigned
    assert g.id == 1
    assert coord.id == 1
    # check we can access coordinates from g
    coord_q = g.coordinates[0]
    assert coord is coord_q
    assert g is coord.genome
    assert g.id == coord_q.genome_id
    # check we get logical behavior on deletion
    sess.delete(coord)
    sess.commit()
    assert len(g.coordinates) == 0
    print(coord.genome)
    sess.delete(g)
    sess.commit()
    with pytest.raises(sqlalchemy.exc.InvalidRequestError):
        sess.add(coord)


def test_coordinate_constraints():
    """Check the coordinate constraints"""
    sess = mk_memory_session()
    g = Genome()

    # should be ok
    coors = Coordinate(length=30, seqid='abc', genome=g)
    coors2 = Coordinate(length=400, seqid='abcd', genome=g)
    sess.add_all([coors, coors2])
    sess.commit()

    # start/end constraints
    coors_bad1 = Coordinate(length=-12, seqid='abc', genome=g)
    coors_bad2 = Coordinate(genome=g)
    with pytest.raises(IntegrityError):
        sess.add(coors_bad1)  # start below 1
        sess.commit()
    sess.rollback()
    with pytest.raises(IntegrityError):
        sess.add(coors_bad2)  # no seqid
        sess.commit()


def test_coordinate_insert():
    """Test what happens when we insert two coordinates"""
    sess = mk_memory_session()
    g = Genome()
    coords = Coordinate(length=10, seqid='abc', genome=g)
    coords2 = Coordinate(length=11, seqid='def', genome=g)
    sl = SuperLocus()
    f0 = Feature(coordinate=coords)
    f1 = Feature(coordinate=coords2)
    # should be ok
    sess.add_all([g, sl, coords, coords2, f0, f1])
    assert f0.coordinate.length == 10
    assert f1.coordinate.length == 11


def test_many2many_with_features():
    """Test the many2many tables association_transcript_piece_to_feature and
    association_protein_to_feature
    """
    sl = SuperLocus()
    # one transcript, multiple proteins
    piece0 = TranscriptPiece()
    slated0 = Protein(super_locus=sl)
    slated1 = Protein(super_locus=sl)
    # features representing alternative start codon for proteins on one transcript
    feat0_tss = Feature(transcript_pieces=[piece0])
    feat2_stop = Feature(proteins=[slated0, slated1])
    feat3_start = Feature(proteins=[slated0])
    # test multi features per protein worked
    assert len(slated0.features) == 2
    # test mutli protein per feature worked
    assert len(feat2_stop.proteins) == 2
    assert len(feat3_start.proteins) == 1
    assert len(feat0_tss.proteins) == 0


def test_feature_has_its_things():
    """Test if properties of Feature table are correct and constraints are enforced"""
    sess = mk_memory_session()
    # test feature with nothing much set
    g = Genome()
    c = Coordinate(length=30, seqid='abc', genome=g)
    f = Feature(coordinate=c,
                start=1,
                end=30,
                start_is_biological_start=True,
                end_is_biological_end=True,
                is_plus_strand=True)
    sess.add_all([f, c])
    sess.commit()

    assert f.source is None
    assert f.score is None
    # test feature with
    f1 = Feature(coordinate=c,
                 start=3,
                 end=-1,
                 start_is_biological_start=True,
                 end_is_biological_end=True,
                 is_plus_strand=False)
    assert not f1.is_plus_strand
    assert f1.start == 3
    assert f1.end == -1
    sess.add(f1)
    sess.commit()

    # test too low of start / end coordinates raise an error
    f_should_fail = Feature(coordinate=c,
                            start=-5,
                            end=10,
                            start_is_biological_start=True,
                            end_is_biological_end=True,
                            is_plus_strand=False)
    sess.add(f_should_fail)
    with pytest.raises(sqlalchemy.exc.IntegrityError):
        sess.commit()
    sess.rollback()
    f_should_fail = Feature(coordinate=c,
                            start=5,
                            end=-2,
                            start_is_biological_start=True,
                            end_is_biological_end=True,
                            is_plus_strand=False)
    sess.add(f_should_fail)
    with pytest.raises(sqlalchemy.exc.IntegrityError):
        sess.commit()
    sess.rollback()

    # test wrong class as parameter
    with pytest.raises(KeyError):
        f2 = Feature(coordinate=f)

    f2 = Feature(coordinate=c,
                 start=1,
                 end=3,
                 start_is_biological_start=True,
                 end_is_biological_end=True,
                 is_plus_strand=-1)  # note that 0, and 1 are accepted
    sess.add(f2)
    with pytest.raises(sqlalchemy.exc.StatementError):
        sess.commit()
    sess.rollback()

    # check 'phase is NULL or (not start_is_biological_start or phase = 0)'
    f = Feature(coordinate=c,
                start=1,
                end=3,
                start_is_biological_start=True,
                end_is_biological_end=True,
                is_plus_strand=True,
                phase=1)
    sess.add(f)
    with pytest.raises(sqlalchemy.exc.IntegrityError):
        sess.commit()
    sess.rollback()


def test_partially_remove_coordinate():
    """Add two coordinates to an annotated genome and test if everything ends
    up valid.
    """
    sess = mk_memory_session()
    g = Genome()
    place_holder = Genome()
    coord0 = Coordinate(length=30, seqid='abc', genome=g)
    coord1 = Coordinate(length=330, seqid='def', genome=g)
    sess.add_all([g, coord0, coord1])
    sess.commit()
    assert len(g.coordinates) == 2
    g.coordinates.remove(coord0)
    coord0.genome = place_holder  # else we'll fail the not NULL constraint
    sess.commit()
    # removed from g
    assert len(g.coordinates) == 1
    # but still in table
    assert len(sess.query(Coordinate).all()) == 2


# section: api
def test_transcript_piece_unique_constraints():
    """Add transcript pieces in valid and invalid configurations and test for
    valid outcomes.
    """
    sess = mk_memory_session()
    sl = SuperLocus()
    transcript0 = Transcript(super_locus=sl)
    transcript1 = Transcript(super_locus=sl)

    # test if same position for different transcript_id goes through
    piece_tr0_pos0 = TranscriptPiece(transcript=transcript0, position=0)
    piece_tr1_pos0 = TranscriptPiece(transcript=transcript1, position=0)
    sess.add_all([transcript0, piece_tr0_pos0, piece_tr1_pos0])
    sess.commit()

    # same transcibed_id but different position
    piece_tr0_pos1 = TranscriptPiece(transcript=transcript0, position=1)
    sess.add(piece_tr0_pos1)
    sess.commit()

    # test if unique constraint works
    piece_tr0_pos1_2nd = TranscriptPiece(transcript=transcript0, position=1)
    sess.add(piece_tr0_pos1_2nd)
    with pytest.raises(IntegrityError):
        sess.commit()


def test_order_pieces():
    """Add transcript pieces that consists of one or many features to the db and test
    if the order for the transcript pieces and the features is returned according to the
    position property instead of db insertion order.
    """
    sess = mk_memory_session()
    g = Genome(species='Athaliana', version='1.2', acquired_from='Phytozome12')
    coor = Coordinate(seqid='a', length=1000, genome=g)
    sess.add_all([g, coor])
    sess.commit()
    # setup one transcript handler with pieces
    sl, sl_h = setup_data_handler(SuperLocusHandlerBase, SuperLocus)
    t, t_h = setup_data_handler(TranscriptHandlerBase, Transcript, super_locus=sl)
    # insert in wrong order
    piece1 = TranscriptPiece(position=1)
    piece0 = TranscriptPiece(position=0)
    piece2 = TranscriptPiece(position=2)
    t.transcript_pieces = [piece0, piece1, piece2]
    sess.add_all([t, piece1, piece0, piece2])
    sess.commit()
    # see if they can be ordered as expected overall
    op = t_h.sorted_pieces
    print([piece0, piece1, piece2], 'expected')
    print(op, 'sorted')
    assert op == [piece0, piece1, piece2]


def test_fasta_import():
    """Import and test coordinate information from fasta files"""

    def import_fasta(path):
        controller = ImportController(database_path='sqlite:///:memory:')
        controller.add_sequences(path)
        return controller

    # test import of multiple sequences from one file
    controller = import_fasta('testdata/basic_sequences.fa')
    coords = controller.session.query(Coordinate).all()
    assert len(coords) == 5
    assert coords[0].seqid == '1'
    assert coords[0].length == 405
    assert coords[0].sha1 == 'dc6f3ba2b0c08f7d08053837b810f86cbaa06f38'
    assert coords[0].sequence == 'N' * 405
    assert coords[1].seqid == 'abc'
    assert coords[1].length == 808
    assert coords[1].sequence == 'AAGGCCTT' * 101
    assert coords[2].seqid == 'test123'
    assert coords[2].length == 100
    assert coords[2].sequence == 'A' * 100


def test_import_multiple_genomes():
    controller = ImportController(database_path='sqlite:///:memory:')
    InsertCounterHolder.sync_counters_with_db(controller.session)
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)
    query = controller.session.query

    n_features_1 = query(Feature).count()
    n_coords_1 = query(Coordinate).count()
    n_genomes_1 = query(Genome).count()
    max_id_1 = query(func.max(Feature.id)).one()[0]

    # add two more genomes
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)

    n_features_3 = query(Feature).count()
    n_coords_3 = query(Coordinate).count()
    n_genomes_3 = query(Genome).count()
    assert n_features_1 * 3 == n_features_3
    assert n_coords_1 * 3 == n_coords_3
    assert n_genomes_1 * 3 == n_genomes_3

    max_id_3 = query(func.max(Feature.id)).one()[0]
    assert max_id_1 * 3 == max_id_3

    # test transcript and super locus relation for a super loci on the plus strand
    # we have 8 super loci in each genome
    transcripts_sl_1 = query(SuperLocus).filter(SuperLocus.id == 1).one().transcripts
    transcripts_sl_2 = query(SuperLocus).filter(SuperLocus.id == 9).one().transcripts
    transcripts_sl_3 = query(SuperLocus).filter(SuperLocus.id == 17).one().transcripts
    assert len(transcripts_sl_1) == len(transcripts_sl_2) == len(transcripts_sl_3)

    # test protein and super locus relation
    proteins_sl_1 = query(SuperLocus).filter(SuperLocus.id == 1).one().proteins
    proteins_sl_2 = query(SuperLocus).filter(SuperLocus.id == 9).one().proteins
    proteins_sl_3 = query(SuperLocus).filter(SuperLocus.id == 17).one().proteins
    assert len(proteins_sl_1) == len(proteins_sl_2) == len(proteins_sl_3)

    # test transcript and super locus relation for a super loci on the minus strand
    transcripts_sl_1 = query(SuperLocus).filter(SuperLocus.id == 6).one().transcripts
    transcripts_sl_2 = query(SuperLocus).filter(SuperLocus.id == 14).one().transcripts
    transcripts_sl_3 = query(SuperLocus).filter(SuperLocus.id == 22).one().transcripts
    assert len(transcripts_sl_1) == len(transcripts_sl_2) == len(transcripts_sl_3)

    # test protein and super locus relation
    proteins_sl_1 = query(SuperLocus).filter(SuperLocus.id == 6).one().proteins
    proteins_sl_2 = query(SuperLocus).filter(SuperLocus.id == 14).one().proteins
    proteins_sl_3 = query(SuperLocus).filter(SuperLocus.id == 22).one().proteins
    assert len(proteins_sl_1) == len(proteins_sl_2) == len(proteins_sl_3)


def test_dummyloci_errors():
    """Tests if all errors generated for dummyloci{.gff|.fa} are correct"""

    def error_in_list(error, error_list):
        """searches for the error in a list. removes the error if found.
        error should be a dict and error list a list of orm objects"""
        for e in error_list[:]:  # make a copy at each iteration so we avoid weird errors
            if (error['coord_id'] == e.coordinate.id and error['is_plus_strand'] == e.is_plus_strand
                    and error['start'] == e.start and error['end'] == e.end
                    and error['type'] == e.type.value):
                error_list.remove(e)
                return True
        return False

    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)
    error_types = [t.value for t in types.Errors]
    errors = controller.session.query(Feature).filter(Feature.type.in_(error_types)).all()
    coords = controller.session.query(Coordinate).all()

    # test case 1 - see gff file for more documentation
    # two identical error bars after cds for aligned exon/cds pair
    error = {
        'coord_id': coords[0].id,
        'is_plus_strand': True,
        'start': 120,
        'end': 499,
        'type': types.MISSING_UTR_3P
    }
    assert error_in_list(error, errors)
    assert error_in_list(error, errors)
    error = {
        'coord_id': coords[0].id,
        'is_plus_strand': True,
        'start': 0,
        'end': 110,
        'type': types.MISSING_UTR_5P
    }
    assert error_in_list(error, errors)
    error = {
        'coord_id': coords[0].id,
        'is_plus_strand': True,
        'start': 0,
        'end': 110,
        'type': types.MISSING_START_CODON
    }
    assert error_in_list(error, errors)

    # test case 2
    error = {
        'coord_id': coords[0].id,
        'is_plus_strand': True,
        'start': 499,
        'end': 1099,
        'type': types.EMPTY_SUPER_LOCUS
    }
    assert error_in_list(error, errors)

    # test case 4 (test case 3 is without errors)
    error = {
        'coord_id': coords[0].id,
        'is_plus_strand': True,
        'start': 1499,
        'end': 1619,
        'type': types.MISSING_START_CODON
    }
    assert error_in_list(error, errors)

    #### Coordinate 1 ####

    # test case 6
    error = {
        'coord_id': coords[1].id,
        'is_plus_strand': True,
        'start': 424,
        'end': 524,
        'type': types.WRONG_PHASE_5P
    }
    assert error_in_list(error, errors)
    error = {
        'coord_id': coords[1].id,
        'is_plus_strand': True,
        'start': 540,
        'end': 579,
        'type': types.TOO_SHORT_INTRON
    }
    assert error_in_list(error, errors)

    # test case 7
    error = {
        'coord_id': coords[1].id,
        'is_plus_strand': False,
        'start': 1449,
        'end': 1349,
        'type': types.MISSING_UTR_5P
    }
    assert error_in_list(error, errors)
    error = {
        'coord_id': coords[1].id,
        'is_plus_strand': False,
        'start': 973,
        'end': -1,
        'type': types.MISMATCHED_PHASE_3P
    }
    assert error_in_list(error, errors)

    # test case 8
    error = {
        'coord_id': coords[1].id,
        'is_plus_strand': False,
        # inclusive start of the minus strand from beginning of the sequence
        'start': 1755,
        'end': 1724,  # exclusive end of intergenic region
        'type': types.MISSING_START_CODON
    }
    assert error_in_list(error, errors)
    error = {
        'coord_id': coords[1].id,
        'is_plus_strand': False,
        'start': 1649,
        'end': 1548,
        'type': types.OVERLAPPING_EXONS
    }
    assert error_in_list(error, errors)

    # test that we don't have any errors we don't expect
    assert not errors


def test_case_1():
    """Confirm the existence of all features of test case 1 of dummyloci.gff except
    for error features, which are tested in test_dummyloci_errors().
    Does not test the exact ids or exactly matching object relationships."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)
    query = controller.session.query

    sl = query(SuperLocus).filter(SuperLocus.given_name == 'gene0').one()
    sl_h = SuperLocusHandlerBase(sl)

    super_locus = SuperLocus(given_name='gene0', type=types.SuperLocusAll.gene)
    assert matching_super_loci(sl, super_locus)

    coords = query(Coordinate).all()

    # confirm exisistence of all objects where things could go wrong
    # not testing for TranscriptPieces as trans-splicing is currently not implemented
    # above db level and one piece has to exist for the sl_h.features query to work
    sl_objects = list(sl_h.features) + sl_h.data.transcripts + sl_h.data.proteins

    # first transcript
    transcript = Transcript(given_name='x1', type=types.TranscriptLevelAll.mRNA, super_locus=sl)
    assert orm_object_in_list(transcript, sl_objects)

    protein = Protein(given_name='x1.p', super_locus=sl)
    assert orm_object_in_list(protein, sl_objects)

    feature = Feature(given_name='x1',
                      type=types.OnSequence.transcript_feature,
                      start=0,
                      end=120,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.coding,
                      start=10,
                      end=120,
                      start_is_biological_start=True,
                      end_is_biological_end=False,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.intron,
                      start=21,
                      end=110,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)

    # second transcript
    transcript = Transcript(given_name='y1', type=types.TranscriptLevelAll.mRNA, super_locus=sl)
    assert orm_object_in_list(transcript, sl_objects)

    protein = Protein(given_name='y1.p', super_locus=sl)
    assert orm_object_in_list(protein, sl_objects)

    feature = Feature(given_name='y1',
                      type=types.OnSequence.transcript_feature,
                      start=0,
                      end=400,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.coding,
                      start=10,
                      end=301,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.intron,
                      start=21,
                      end=110,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.intron,
                      start=120,
                      end=200,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)

    # third transcript
    transcript = Transcript(given_name='z1', type=types.TranscriptLevelAll.mRNA, super_locus=sl)
    assert orm_object_in_list(transcript, sl_objects)

    protein = Protein(given_name='z1.p', super_locus=sl)
    assert orm_object_in_list(protein, sl_objects)

    feature = Feature(given_name='z1',
                      type=types.OnSequence.transcript_feature,
                      start=110,
                      end=120,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.coding,
                      start=110,
                      end=120,
                      start_is_biological_start=False,
                      end_is_biological_end=False,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)

    # test if we have no extra objects
    assert not sl_objects


def test_case_8():
    """Analogous to test_case_1()"""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)
    query = controller.session.query

    sl = query(SuperLocus).\
        filter(SuperLocus.given_name == 'gene_overlapping_exons_missing_start').one()
    sl_h = SuperLocusHandlerBase(sl)

    super_locus = SuperLocus(given_name='gene_overlapping_exons_missing_start',
                             type=types.SuperLocusAll.gene)
    assert matching_super_loci(sl, super_locus)

    coords = query(Coordinate).all()

    sl_objects = list(sl_h.features) + sl_h.data.transcripts + sl_h.data.proteins

    # first transcript
    transcript = Transcript(given_name='x8', type=types.TranscriptLevelAll.mRNA, super_locus=sl)
    assert orm_object_in_list(transcript, sl_objects)

    protein = Protein(given_name='x8.p', super_locus=sl)
    assert orm_object_in_list(protein, sl_objects)

    feature = Feature(given_name='x8',
                      type=types.OnSequence.transcript_feature,
                      start=1749,
                      end=1548,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0,
                      coordinate=coords[1])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.coding,
                      start=1724,
                      end=1573,
                      start_is_biological_start=False,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0,
                      coordinate=coords[1])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.intron,
                      start=1718,
                      end=1649,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0,
                      coordinate=coords[1])
    assert orm_object_in_list(feature, sl_objects)

    # second transcript
    transcript = Transcript(given_name='y8', type=types.TranscriptLevelAll.mRNA, super_locus=sl)
    assert orm_object_in_list(transcript, sl_objects)

    protein = Protein(given_name='y8.p', super_locus=sl)
    assert orm_object_in_list(protein, sl_objects)

    feature = Feature(given_name='y8',
                      type=types.OnSequence.transcript_feature,
                      start=1749,
                      end=1548,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0,
                      coordinate=coords[1])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.coding,
                      start=1729,
                      end=1699,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0,
                      coordinate=coords[1])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.intron,
                      start=1678,
                      end=1599,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0,
                      coordinate=coords[1])
    assert orm_object_in_list(feature, sl_objects)

    # test if we have no extra objects
    assert not sl_objects


def test_gff_gen():
    gff_organizer = OrganizedGFFEntries('testdata/testerSl.gff3')
    x = list(gff_organizer._gff_gen())
    assert len(x) == 103
    assert x[0].type == 'region'
    assert x[-1].type == 'CDS'


def test_gff_useful_gen():
    gff_organizer = OrganizedGFFEntries('testdata/testerSl.gff3')
    x = list(gff_organizer._useful_gff_entries())
    assert len(x) == 100  # should drop the region entry
    assert x[0].type == 'gene'
    assert x[-1].type == 'CDS'


def test_gff_grouper():
    gff_organizer = OrganizedGFFEntries('testdata/testerSl.gff3')
    gff_organizer.load_organized_entries()
    n_genes_seqid = {'NC_015438.2': 2, 'NC_015439.2': 2, 'NC_015440.2': 1}
    for seqid, count in n_genes_seqid.items():
        assert len(gff_organizer.organized_entries[seqid]) == count
        for group in gff_organizer.organized_entries[seqid]:
            assert group[0].type == 'gene'


# section: types
def test_enum_non_inheritance():
    allknown = [x.name for x in list(types.AllKnown)]
    allnice = [x.name for x in list(types.AllKeepable)]
    # check that some random bits made it in to all
    assert 'missing_utr_3p' in allknown
    assert 'region' in allknown

    # check that some annoying bits are not in nice set
    for not_nice in ['transcript', 'primary_transcript', 'exon', 'five_prime_UTR', 'CDS']:
        assert not_nice not in allnice
        assert not_nice in allknown

    # check nothing is there twice
    assert len(set(allknown)) == len(allknown)


def test_enums_name_val_match():
    for x in types.AllKnown:
        assert x.name == x.value


# section: helpers
def test_key_matching():
    # identical
    known = {'a', 'b', 'c'}
    mapper, is_forward = helpers.two_way_key_match(known, known)
    assert is_forward
    assert mapper('a') == 'a'
    assert isinstance(mapper, helpers.CheckMapper)
    with pytest.raises(KeyError):
        mapper('d')

    # subset, should behave as before
    mapper, is_forward = helpers.two_way_key_match(known, {'c'})
    assert is_forward
    assert mapper('a') == 'a'
    assert isinstance(mapper, helpers.CheckMapper)
    with pytest.raises(KeyError):
        mapper('d')

    # superset, should flip ordering
    mapper, is_forward = helpers.two_way_key_match({'c'}, known)
    assert not is_forward
    assert mapper('a') == 'a'
    assert isinstance(mapper, helpers.CheckMapper)
    with pytest.raises(KeyError):
        mapper('d')

    # other is abbreviated from known
    set1 = {'a.seq', 'b.seq', 'c.seq'}
    set2 = {'a', 'b', 'c'}
    mapper, is_forward = helpers.two_way_key_match(set1, set2)
    assert is_forward
    assert mapper('a') == 'a.seq'
    assert isinstance(mapper, helpers.DictMapper)
    with pytest.raises(KeyError):
        mapper('d')

    # cannot be safely differentiated
    set1 = {'ab.seq', 'ba.seq', 'c.seq', 'a.seq', 'b.seq'}
    set2 = {'a', 'b', 'c'}
    with pytest.raises(helpers.NonMatchableIDs):
        mapper, is_forward = helpers.two_way_key_match(set1, set2)
        print(mapper.key_vals)


def test_gff_to_seqids():
    x = helpers.get_seqids_from_gff('testdata/testerSl.gff3')
    assert x == {'NC_015438.2', 'NC_015439.2', 'NC_015440.2'}
