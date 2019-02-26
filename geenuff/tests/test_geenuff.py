from dustdas import gffhelper
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import sqlalchemy
import pytest
import os

#import sys
#sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))

from .. import orm
from .. import api
from ..applications import gffimporter
from .. import types
from .. import helpers

# section: annotations_orm
def mk_session(db_path='sqlite:///:memory:'):
    engine = create_engine(db_path, echo=False)
    orm.Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    return Session()


def test_annogenome2sequence_infos_relation():
    sess = mk_session()
    ag = orm.AnnotatedGenome(species='Athaliana', version='1.2', acquired_from='Phytozome12')
    sequence_info = orm.SequenceInfo(annotated_genome=ag)
    assert ag is sequence_info.annotated_genome
    # actually put everything in db
    sess.add(sequence_info)
    sess.commit()
    # check primary keys were assigned
    assert ag.id == 1
    assert sequence_info.id == 1
    # check we can access sequence_info from ag
    sequence_info_q = ag.sequence_infos[0]
    assert sequence_info is sequence_info_q
    assert ag is sequence_info.annotated_genome
    assert ag.id == sequence_info_q.annotated_genome_id
    # check we get logical behavior on deletion
    sess.delete(sequence_info)
    sess.commit()
    assert len(ag.sequence_infos) == 0
    print(sequence_info.annotated_genome)
    sess.delete(ag)
    sess.commit()
    with pytest.raises(sqlalchemy.exc.InvalidRequestError):
        sess.add(sequence_info)


def test_coordinate_constraints():
    sess = mk_session()
    coors = orm.Coordinates(start=0, end=30, seqid='abc')
    coors2 = orm.Coordinates(start=0, end=1, seqid='abc')
    coors_bad1 = orm.Coordinates(start=-12, end=30, seqid='abc')
    coors_bad2 = orm.Coordinates(start=100, end=30, seqid='abc')
    coors_bad3 = orm.Coordinates(start=1, end=30)
    # should be ok
    sess.add_all([coors, coors2])
    sess.commit()
    # should cause trouble
    with pytest.raises(sqlalchemy.exc.IntegrityError):
        sess.add(coors_bad1)  # start below 1
        sess.commit()
    sess.rollback()
    with pytest.raises(sqlalchemy.exc.IntegrityError):
        sess.add(coors_bad2)  # end below start
        sess.commit()
    sess.rollback()
    with pytest.raises(sqlalchemy.exc.IntegrityError):
        sess.add(coors_bad3)
        sess.commit()


def test_coordinate_seqinfo_query():
    sess = mk_session()
    ag = orm.AnnotatedGenome()
    si_model = orm.SequenceInfo(annotated_genome=ag)
    slic = api.SequenceInfoHandler()
    slic.add_data(si_model)
    coors = orm.Coordinates(start=1, end=30, seqid='abc', sequence_info=si_model)
    coors2 = orm.Coordinates(start=11, end=330, seqid='def', sequence_info=si_model)
    sl = orm.SuperLocus()
    f0 = orm.Feature(super_locus=sl, coordinates=coors)
    f1 = orm.Feature(super_locus=sl, coordinates=coors2)
    # should be ok
    sess.add_all([coors, coors2, sl, f0, f1])
    assert f0.coordinates.start == 1
    assert f1.coordinates.end == 330
    #assert seq_info is slic.seq_info


def test_many2many_scribed2slated():
    # test transcript to multi proteins
    sl = orm.SuperLocus()
    scribed0 = orm.Transcribed(super_locus=sl)
    slated0 = orm.Translated(super_locus=sl, transcribeds=[scribed0])
    slated1 = orm.Translated(super_locus=sl, transcribeds=[scribed0])
    # test scribed 2 scribed_piece works
    piece0 = orm.TranscribedPiece(transcribed=scribed0)
    assert piece0 == slated0.transcribeds[0].transcribed_pieces[0]
    # test protein to multi transcripts
    assert set(scribed0.translateds) == {slated0, slated1}
    slated2 = orm.Translated(super_locus=sl)
    scribed1 = orm.Transcribed(super_locus=sl)
    scribed2 = orm.Transcribed(super_locus=sl)
    slated2.transcribeds = [scribed1, scribed2]
    assert set(slated2.transcribeds) == {scribed1, scribed2}
    scribed3 = orm.Transcribed(super_locus=sl)
    slated2.transcribeds.append(scribed3)
    assert set(slated2.transcribeds) == {scribed1, scribed2, scribed3}


def test_many2many_with_features():
    sl = orm.SuperLocus()
    # one transcript, multiple proteins
    piece0 = orm.TranscribedPiece(super_locus=sl)
    scribed0 = orm.Transcribed(super_locus=sl, transcribed_pieces=[piece0])
    slated0 = orm.Translated(super_locus=sl, transcribeds=[scribed0])
    slated1 = orm.Translated(super_locus=sl, transcribeds=[scribed0])
    # features representing alternative start codon for proteins on one transcript
    feat0_tss = orm.Feature(super_locus=sl, transcribed_pieces=[piece0])
    feat1_tss = orm.Feature(super_locus=sl, transcribed_pieces=[piece0])
    feat2_stop = orm.Feature(super_locus=sl, translateds=[slated0, slated1])
    feat3_start = orm.Feature(super_locus=sl, translateds=[slated0])
    feat4_start = orm.Feature(super_locus=sl, translateds=[slated1])
    # test they all made it to super locus
    assert len(sl.features) == 5
    # test multi features per translated worked
    assert len(slated0.features) == 2
    # test mutli translated per feature worked
    assert len(feat2_stop.translateds) == 2
    assert len(feat3_start.translateds) == 1
    assert len(feat0_tss.translateds) == 0
    # test we can get to all of this from transcribed
    indirect_features = set()
    for slated in scribed0.translateds:
        for f in slated.features:
            indirect_features.add(f)
    assert len(indirect_features) == 3


def test_feature_has_its_things():
    sess = mk_session()
    # should be ok
    sl = orm.SuperLocus()
    # test feature with nothing much set
    f = orm.Feature(super_locus=sl)
    sess.add(f)
    sess.commit()

    assert f.is_plus_strand is None
    assert f.source is None
    assert f.coordinates is None
    assert f.score is None
    # test feature with
    f1 = orm.Feature(super_locus=sl, is_plus_strand=False, position=3)
    assert not f1.is_plus_strand
    assert f1.position == 3
    # test bad input
    with pytest.raises(KeyError):
        f2 = orm.Feature(super_locus=f)

    f2 = orm.Feature(is_plus_strand=-1)  # note that 0, and 1 are accepted
    sess.add(f2)
    with pytest.raises(sqlalchemy.exc.StatementError):
        sess.commit()
    sess.rollback()


def test_feature_streamlinks():
    sess = mk_session()
    f = orm.Feature(position=1)
    pair = orm.UpDownPair()
    sfA0 = orm.UpstreamFeature(position=2, pairs=[pair])
    sfA1 = orm.DownstreamFeature(position=3, pairs=[pair])
    sess.add_all([f, sfA0, sfA1, pair])
    sess.commit()
    sfA0back = sess.query(orm.UpstreamFeature).first()
    assert sfA0back is sfA0
    sfA1back = sfA0back.pairs[0].downstream
    assert sfA1 is sfA1back
    assert sfA1.pairs[0].upstream is sfA0
    sf_friendless = orm.DownstreamFeature(position=4)
    sess.add(sf_friendless)
    sess.commit()
    downstreams = sess.query(orm.DownstreamFeature).all()
    pre_downlinked = sess.query(orm.UpDownPair).filter(
        orm.DownstreamFeature != None  # todo, isn't there a 'right' way to do this?
    ).all()
    downlinked = [x.downstream for x in pre_downlinked]
    print([(x.position, x.pairs) for x in downstreams])
    print([(x.position, x.pairs[0].upstream) for x in downlinked])
    assert len(downstreams) == 2
    assert len(downlinked) == 1


def test_linking_via_fkey():
    sess = mk_session()
    sfA0 = orm.UpstreamFeature(position=2)
    sfA1 = orm.DownstreamFeature(position=3)
    sess.add_all([sfA0, sfA1])
    sess.commit()
    pair = orm.UpDownPair(upstream_id=sfA0.id, downstream_id=sfA1.id)
    sess.add_all([pair])
    sess.commit()
    assert sfA1.pairs[0].upstream is sfA0
    assert sfA0.pairs[0].downstream is sfA1


def test_delinking_from_oneside():
    sess = mk_session()
    ag = orm.AnnotatedGenome()
    place_holder = orm.AnnotatedGenome()
    si0 = orm.SequenceInfo(annotated_genome=ag)
    si1 = orm.SequenceInfo(annotated_genome=ag)
    sess.add_all([ag, si0, si1])
    sess.commit()
    assert len(ag.sequence_infos) == 2
    ag.sequence_infos.remove(si0)
    si0.annotated_genome = place_holder  # else we'll fail the not NULL constraint
    sess.commit()
    # removed from ag
    assert len(ag.sequence_infos) == 1
    # but still in table
    assert len(sess.query(orm.SequenceInfo).all()) == 2


# section: annotations
def test_copy_over_attr():
    sess = mk_session()
    data_ag, dummy_ag = setup_data_handler(api.AnnotatedGenomeHandler, orm.AnnotatedGenome,
                                           species='mammoth', version='v1.0.3', acquired_from='nowhere')
    odata_ag, other_ag = setup_data_handler(api.AnnotatedGenomeHandler, orm.AnnotatedGenome)

    sess.add_all([data_ag, other_ag.data])
    # make sure copy only copies what it says and nothing else
    dummy_ag.copy_data_attr_to_other(other_ag, copy_only='species')
    assert other_ag.get_data_attribute('species') == 'mammoth'
    assert other_ag.get_data_attribute('version') is None
    sess.add_all([odata_ag, data_ag])
    # make sure do_not_copy excludes what is says and copies the rest
    dummy_ag.copy_data_attr_to_other(other_ag, do_not_copy='acquired_from')
    assert other_ag.get_data_attribute('version') == 'v1.0.3'
    assert other_ag.get_data_attribute('acquired_from') is None
    # make sure everything is copied
    dummy_ag.copy_data_attr_to_other(other_ag)
    assert other_ag.get_data_attribute('acquired_from') == 'nowhere'
    # make sure commit/actual entry works
    # sess.add_all([data_ag, other_ag.data])
    sess.commit()
    assert dummy_ag.get_data_attribute('species') == 'mammoth'
    assert other_ag.get_data_attribute('species') == 'mammoth'
    assert other_ag.get_data_attribute('acquired_from') == 'nowhere'
    assert other_ag.data.id is not None


def test_swap_link_annogenome2seqinfo():
    sess = mk_session()

    ag, agh = setup_data_handler(api.AnnotatedGenomeHandler, orm.AnnotatedGenome)
    ag2, ag2h = setup_data_handler(api.AnnotatedGenomeHandler, orm.AnnotatedGenome)

    si, sih = setup_data_handler(api.SequenceInfoHandler, orm.SequenceInfo, annotated_genome=ag)

    sess.add_all([ag, ag2, si])
    sess.commit()
    assert agh.data.sequence_infos == [sih.data]
    agh.de_link(sih)
    ag2h.link_to(sih)
    sess.commit()
    assert agh.data.sequence_infos == []
    assert ag2h.data.sequence_infos == [sih.data]
    # swap back from sequence info interface
    sih.de_link(ag2h)
    sih.link_to(agh)
    assert agh.data.sequence_infos == [sih.data]
    assert ag2h.data.sequence_infos == []


def test_swap_links_superlocus2ttfs():
    sess = mk_session()

    slc, slch = setup_data_handler(api.SuperLocusHandler, orm.SuperLocus)

    slc2, slc2h = setup_data_handler(api.SuperLocusHandler, orm.SuperLocus)

    scribed, scribedh = setup_data_handler(api.TranscribedHandler, orm.Transcribed,
                                           super_locus=slc)
    slated, slatedh = setup_data_handler(api.TranslatedHandler, orm.Translated,
                                         super_locus=slc, transcribeds=[scribed])
    feature, featureh = setup_data_handler(api.FeatureHandler, orm.Feature,
                                           super_locus=slc)

    sess.add_all([slc, slc2, scribed, slated])
    sess.commit()
    # swapping super locus
    slch.de_link(slatedh)
    slch.de_link(scribedh)
    slch.de_link(featureh)
    slc2h.link_to(slatedh)
    slc2h.link_to(scribedh)
    slc2h.link_to(featureh)
    assert scribed.super_locus is slc2
    assert slated.super_locus is slc2
    assert feature.super_locus is slc2
    assert slc.translateds == []
    assert slc.transcribeds == []
    assert slc.features == []
    # swapping back from transcribed, translated, feature side
    scribedh.de_link(slc2h)
    slatedh.de_link(slc2h)
    featureh.de_link(slc2h)
    scribedh.link_to(slch)
    slatedh.link_to(slch)
    featureh.link_to(slch)
    assert slated.super_locus is slc
    assert scribed.super_locus is slc
    assert feature.super_locus is slc
    assert slc2.translateds == []
    assert slc2.transcribeds == []
    assert slc2.features == []
    sess.commit()


def test_swap_links_t2t2f():
    sess = mk_session()

    slc = orm.SuperLocus()

    scribedpiece, scribedpieceh = setup_data_handler(api.TranscribedPieceHandler,
                                                     orm.TranscribedPiece, super_locus=slc)
    scribed, scribedh = setup_data_handler(api.TranscribedHandler, orm.Transcribed,
                                           super_locus=slc, transcribed_pieces=[scribedpiece])
    slated, slatedh = setup_data_handler(api.TranslatedHandler, orm.Translated, super_locus=slc,
                                         transcribeds=[scribed])
    feature, featureh = setup_data_handler(api.FeatureHandler, orm.Feature, super_locus=slc,
                                           transcribed_pieces=[scribedpiece])

    sess.add_all([slc, scribed, scribedpiece, slated, feature])
    sess.commit()

    assert scribed.translateds == [slated]
    assert slated.transcribeds == [scribed]
    assert scribedpiece.transcribed == scribed
    assert feature.transcribed_pieces == [scribedpiece]
    # de_link / link_to from scribed side
    scribedh.de_link(slatedh)
    scribedpieceh.de_link(featureh)
    assert slated.transcribeds == []
    assert scribed.translateds == []
    assert feature.transcribed_pieces == []

    scribedh.link_to(slatedh)
    assert scribed.translateds == [slated]
    assert slated.transcribeds == [scribed]
    # de_link / link_to from slated side
    slatedh.de_link(scribedh)
    assert slated.transcribeds == []
    assert scribed.translateds == []
    slatedh.link_to(scribedh)
    slatedh.link_to(featureh)
    assert scribed.translateds == [slated]
    assert slated.transcribeds == [scribed]
    assert feature.translateds == [slated]
    # mod links from feature side
    featureh.de_link(slatedh)
    featureh.link_to(scribedpieceh)
    assert slated.features == []
    assert scribed.transcribed_pieces[0].features == [feature]
    sess.commit()


def test_updownhandler_links():
    sess = mk_session()
    coor_old = orm.Coordinates(start=1, end=1000, seqid='a')
    coor_new = orm.Coordinates(start=1, end=100, seqid='a')
    slc = orm.SuperLocus()
    scribedpiece, scribedpieceh = setup_data_handler(api.TranscribedPieceHandler,
                                                     orm.TranscribedPiece, super_locus=slc)
    scribed, scribedh = setup_data_handler(api.TranscribedHandler, orm.Transcribed,
                                           super_locus=slc, transcribed_pieces=[scribedpiece])
    up, uph = setup_data_handler(api.UpstreamFeatureHandler, orm.UpstreamFeature, super_locus=slc,
                                 transcribed_pieces=[scribedpiece], coordinates=coor_old)
    up2, up2h = setup_data_handler(api.UpstreamFeatureHandler, orm.UpstreamFeature, super_locus=slc,
                                   transcribed_pieces=[scribedpiece], coordinates=coor_new)

    down, downh = setup_data_handler(api.DownstreamFeatureHandler, orm.DownstreamFeature,
                                     coordinates=coor_old)

    pair, pairh = setup_data_handler(api.UpDownPairHandler, orm.UpDownPair, transcribed=scribed,
                                     upstream=up, downstream=down)
    sess.add_all([up, up2, down, slc, coor_old, coor_new, pair])
    sess.commit()
    assert up2.pairs == []
    assert up.pairs[0].downstream == down
    pairh.de_link(uph)
    pairh.link_to(up2h)
    assert up2.pairs[0].downstream == down
    assert up.pairs == []


def setup_data_handler(handler_type, data_type, **kwargs):
    data = data_type(**kwargs)
    handler = handler_type()
    handler.add_data(data)
    return data, handler


def test_replacelinks():
    sess = mk_session()
    slc = orm.SuperLocus()

    scribedpiece, scribedpieceh = setup_data_handler(api.TranscribedPieceHandler,
                                                     orm.TranscribedPiece, super_locus=slc)
    assert scribedpiece.super_locus is slc
    slated, slatedh = setup_data_handler(api.TranslatedHandler, orm.Translated, super_locus=slc)
    f0, f0h = setup_data_handler(api.FeatureHandler, orm.Feature, super_locus=slc,
                                 translateds=[slated])

    f1, f1h = setup_data_handler(api.FeatureHandler, orm.Feature, super_locus=slc,
                                 translateds=[slated])

    f2, f2h = setup_data_handler(api.FeatureHandler, orm.Feature, super_locus=slc,
                                 translateds=[slated])

    sess.add_all([slc, scribedpiece, slated, f0, f1, f2])
    sess.commit()
    assert len(slated.features) == 3
    slatedh.replace_selflinks_w_replacementlinks(replacement=scribedpieceh, to_replace=['features'])

    assert len(slated.features) == 0
    assert len(scribedpiece.features) == 3


def test_order_pieces():
    sess = mk_session()
    ag = orm.AnnotatedGenome(species='Athaliana', version='1.2', acquired_from='Phytozome12')
    sequence_info = orm.SequenceInfo(annotated_genome=ag)
    coor = orm.Coordinates(seqid='a', start=1, end=1000, sequence_info=sequence_info)
    sess.add_all([ag, sequence_info, coor])
    sess.commit()
    # setup one transcribed handler with pieces
    scribed, scribedh = setup_data_handler(api.TranscribedHandler, orm.Transcribed)
    ti = api.TranscriptInterpBase(transcript=scribedh, session=sess)
    piece1 = orm.TranscribedPiece()
    piece0 = orm.TranscribedPiece()
    piece2 = orm.TranscribedPiece()
    scribed.transcribed_pieces = [piece0, piece1, piece2]
    sess.add_all([scribed, piece0, piece1, piece2])
    sess.commit()
    # setup some paired features
    feature0u = orm.UpstreamFeature(transcribed_pieces=[piece0], coordinates=coor, position=100, given_id='0u',
                                    is_plus_strand=True)
    feature1d = orm.DownstreamFeature(transcribed_pieces=[piece1], coordinates=coor, position=1, given_id='1d',
                                      is_plus_strand=True)
    feature1u = orm.UpstreamFeature(transcribed_pieces=[piece1], coordinates=coor, position=100, given_id='1u',
                                    is_plus_strand=True)
    feature2d = orm.DownstreamFeature(transcribed_pieces=[piece2], coordinates=coor, position=1, given_id='2d',
                                      is_plus_strand=True)
    pair01 = orm.UpDownPair(upstream=feature0u, downstream=feature1d, transcribed=scribed)
    pair12 = orm.UpDownPair(upstream=feature1u, downstream=feature2d, transcribed=scribed)
    # check getting upstream link
    upstream_link = ti.get_upstream_link(piece1)
    assert upstream_link is pair01
    upstream = upstream_link.upstream
    assert upstream is feature0u
    assert upstream.transcribed_pieces == [piece0]
    # check getting downstream link
    downstream_link = ti.get_downstream_link(piece1)
    assert downstream_link is pair12
    downstream = downstream_link.downstream
    assert downstream is feature2d
    assert downstream.transcribed_pieces == [piece2]
    # and see if they can be ordered as expected overall
    op = ti.sort_pieces()
    print([piece0, piece1, piece2], 'expected')
    print(op, 'sorted')
    assert op == [piece0, piece1, piece2]
    # and finally, order features by piece
    fully_sorted = ti.sort_all()
    expected = [[feature0u],
                [feature1d, feature1u],
                [feature2d]]
    assert fully_sorted == expected
    # finally make it circular, and make sure it throws an error
    feature2u = orm.UpstreamFeature(transcribed_pieces=[piece2], coordinates=coor, position=100, given_id='2u',
                                    is_plus_strand=True)
    feature0d = orm.DownstreamFeature(transcribed_pieces=[piece0], coordinates=coor, position=1, given_id='0d',
                                      is_plus_strand=True)
    pair20 = orm.UpDownPair(upstream=feature2u, downstream=feature0d, transcribed=scribed)
    sess.add(pair20)
    sess.commit()
    with pytest.raises(api.IndecipherableLinkageError):
        ti.sort_pieces()


class TransspliceDemoData(object):
    def __init__(self, sess):
        # setup two transitions:
        # 1) scribed - [->[TSS(A),START(B),TDSS(C{->F}),TTS(D)], ->[TSS(E), <<slice>>> TASS(F),STOP(G),TTS(H)]]
        # 2) scribedflip - [->[TSS(A),START(B),TDSS(C{->F'}),TTS(D)], <-[TTS(H'), <<slice>> STOP(G'),TASS(F'),TSS(E')]]
        self.old_coor = orm.Coordinates(seqid='a', start=1, end=2000)
        self.sl, self.slh = setup_data_handler(api.SuperLocusHandler, orm.SuperLocus)
        self.scribed, self.scribedh = setup_data_handler(api.TranscribedHandler, orm.Transcribed,
                                                         super_locus=self.sl)
        self.scribedflip, self.scribedfliph = setup_data_handler(api.TranscribedHandler,
                                                                 orm.Transcribed,
                                                                 super_locus=self.sl)

        self.ti = api.TranscriptInterpBase(transcript=self.scribedh, session=sess)
        self.tiflip = api.TranscriptInterpBase(transcript=self.scribedfliph, session=sess)

        self.pieceA2D = orm.TranscribedPiece(super_locus=self.sl)
        self.pieceA2Dp = orm.TranscribedPiece(super_locus=self.sl)
        self.pieceE2H = orm.TranscribedPiece(super_locus=self.sl)
        self.pieceEp2Hp = orm.TranscribedPiece(super_locus=self.sl)
        self.scribed.transcribed_pieces = [self.pieceA2D, self.pieceE2H]
        self.scribedflip.transcribed_pieces = [self.pieceA2Dp, self.pieceEp2Hp]
        # pieceA2D features

        self.fA = orm.Feature(coordinates=self.old_coor, position=10, given_id='A',
                              is_plus_strand=True, super_locus=self.sl,
                              type=types.TRANSCRIBED, bearing=types.START)
        self.fB = orm.Feature(coordinates=self.old_coor, position=20, given_id='B',
                              is_plus_strand=True, super_locus=self.sl, type=types.CODING,
                              bearing=types.START)

        self.fC = orm.Feature(coordinates=self.old_coor, position=30, given_id='C',
                              is_plus_strand=True, super_locus=self.sl,
                              type=types.TRANS_INTRON, bearing=types.START)
        self.fD = orm.Feature(coordinates=self.old_coor, position=40, given_id='D',
                              is_plus_strand=True, super_locus=self.sl,
                              type=types.TRANSCRIBED, bearing=types.END)
        self.fADs0 = orm.UpstreamFeature(coordinates=self.old_coor, position=40, given_id='ADs0',
                                         is_plus_strand=True, super_locus=self.sl,
                                         type=types.TRANS_INTRON, bearing=types.CLOSE_STATUS)
        self.fADs1 = orm.UpstreamFeature(coordinates=self.old_coor, position=40, given_id='ADs1',
                                         is_plus_strand=True, super_locus=self.sl,
                                         type=types.CODING, bearing=types.CLOSE_STATUS)
        # pieceE2H features
        self.fEHs0 = orm.DownstreamFeature(coordinates=self.old_coor, position=910, given_id='EHs0',
                                           is_plus_strand=True, super_locus=self.sl,
                                           type=types.TRANS_INTRON, bearing=types.OPEN_STATUS)
        self.fEHs1 = orm.DownstreamFeature(coordinates=self.old_coor, position=910, given_id='EHs1',
                                           is_plus_strand=True, super_locus=self.sl,
                                           type=types.CODING, bearing=types.OPEN_STATUS)
        self.fE = orm.Feature(coordinates=self.old_coor, position=910, given_id='E',
                              is_plus_strand=True, super_locus=self.sl,
                              type=types.TRANSCRIBED, bearing=types.START)
        self.fF = orm.Feature(coordinates=self.old_coor, position=920, given_id='F',
                              super_locus=self.sl, is_plus_strand=True,
                              type=types.TRANS_INTRON, bearing=types.END)
        self.fG = orm.Feature(coordinates=self.old_coor, position=930, given_id='G',
                              is_plus_strand=True, super_locus=self.sl, type=types.CODING,
                              bearing=types.END)
        self.fH = orm.Feature(coordinates=self.old_coor, position=940, given_id='H',
                              is_plus_strand=True, super_locus=self.sl,
                              type=types.TRANSCRIBED, bearing=types.END)
        # pieceEp2Hp features
        self.fEHps0 = orm.DownstreamFeature(coordinates=self.old_coor, position=940, given_id='EHsp0',
                                            is_plus_strand=False, super_locus=self.sl,
                                            type=types.TRANS_INTRON, bearing=types.OPEN_STATUS)
        self.fEHps1 = orm.DownstreamFeature(coordinates=self.old_coor, position=940, given_id='EHsp1',
                                            is_plus_strand=False, super_locus=self.sl,
                                            type=types.CODING, bearing=types.OPEN_STATUS)
        self.fEp = orm.Feature(coordinates=self.old_coor, position=940, given_id='Ep',
                               is_plus_strand=False, super_locus=self.sl,
                               type=types.TRANSCRIBED, bearing=types.START)
        self.fFp = orm.Feature(coordinates=self.old_coor, position=930, given_id='Fp',
                               super_locus=self.sl, bearing=types.END,
                               is_plus_strand=False, type=types.TRANS_INTRON)
        self.fGp = orm.Feature(coordinates=self.old_coor, position=920, given_id='Gp',
                               is_plus_strand=False, super_locus=self.sl, type=types.CODING,
                               bearing=types.END)
        self.fHp = orm.Feature(coordinates=self.old_coor, position=910, given_id='Hp',
                               is_plus_strand=False, super_locus=self.sl,
                               type=types.TRANSCRIBED, bearing=types.END)

        self.pieceA2D.features = [self.fA, self.fB, self.fC, self.fD, self.fADs0, self.fADs1]
        self.pieceA2Dp.features = [self.fA, self.fB, self.fC, self.fD, self.fADs0, self.fADs1]
        self.pieceE2H.features = [self.fE, self.fF, self.fG, self.fH, self.fEHs0, self.fEHs1]
        self.pieceEp2Hp.features = [self.fEp, self.fFp, self.fGp, self.fHp, self.fEHps0, self.fEHps1]
        self.pairADEH0 = orm.UpDownPair(upstream=self.fADs0, downstream=self.fEHs0,
                                        transcribed=self.scribed)
        self.pairADEH1 = orm.UpDownPair(upstream=self.fADs1, downstream=self.fEHs1,
                                        transcribed=self.scribed)
        self.pairADEHp0 = orm.UpDownPair(upstream=self.fADs0, downstream=self.fEHps0,
                                         transcribed=self.scribedflip)
        self.pairADEHp1 = orm.UpDownPair(upstream=self.fADs1, downstream=self.fEHps1,
                                         transcribed=self.scribedflip)
        sess.add_all([self.sl, self.pairADEH0, self.pairADEH1, self.pairADEHp0, self.pairADEHp1])
        sess.commit()

    def make_all_handlers(self):
        self.slh.make_all_handlers()


def test_transition_transsplice():
    sess = mk_session()
    d = TransspliceDemoData(sess)  # setup _d_ata
    d.make_all_handlers()
    # forward pass, same sequence, two pieces
    ti_transitions = list(d.ti.transition_5p_to_3p())
    # from transition gen: 0 -> aligned_Features, 1 -> status copy
    assert [set(x[0]) for x in ti_transitions] == [{d.fA}, {d.fB}, {d.fC}, {d.fD}, {d.fADs1, d.fADs0},
                                                   {d.fEHs0, d.fEHs1}, {d.fE}, {d.fF}, {d.fG}, {d.fH}]
    print([x[1].genic for x in ti_transitions])
    assert [x[1].genic for x in ti_transitions] == [bool(x) for x in [1, 1, 1, 0, 0, 0, 1, 1, 1, 0]]
    assert [x[1].in_translated_region for x in ti_transitions] == [bool(x) for x in [0, 1, 1, 1, 0, 1, 1, 1, 0, 0]]
    assert [x[1].in_trans_intron for x in ti_transitions] == [bool(x) for x in [0, 0, 1, 1, 0, 1, 1, 0, 0, 0]]
    # forward, then backward pass, same sequence, two pieces
    ti_transitions = list(d.tiflip.transition_5p_to_3p())
    assert [set(x[0]) for x in ti_transitions] == [{d.fA}, {d.fB}, {d.fC}, {d.fD}, {d.fADs1, d.fADs0},
                                                   {d.fEHps0, d.fEHps1}, {d.fEp}, {d.fFp}, {d.fGp}, {d.fHp}]
    print([x[1].genic for x in ti_transitions])
    assert [x[1].genic for x in ti_transitions] == [bool(x) for x in [1, 1, 1, 0, 0, 0, 1, 1, 1, 0]]
    assert [x[1].in_translated_region for x in ti_transitions] == [bool(x) for x in [0, 1, 1, 1, 0, 1, 1, 1, 0, 0]]
    assert [x[1].in_trans_intron for x in ti_transitions] == [bool(x) for x in [0, 0, 1, 1, 0, 1, 1, 0, 0, 0]]


# section: gff_2_annotations
def test_data_frm_gffentry():
    #sess = mk_session()
    controller = gffimporter.ImportControl(database_path='sqlite:///:memory:', err_path=None)

    sess = controller.session
    ag = orm.AnnotatedGenome()
    slice, sliceh = setup_data_handler(gffimporter.SequenceInfoHandler, orm.SequenceInfo,
                                       annotated_genome=ag)
    coors = orm.Coordinates(seqid='NC_015438.2', start=1, end=100000, sequence_info=slice)
    sess.add_all([ag, slice, coors])
    sess.commit()
    sliceh.mk_mapper()  # todo, why doesn't this work WAS HERE
    gene_string = 'NC_015438.2\tGnomon\tgene\t4343\t5685\t.\t+\t.\tID=gene0;Dbxref=GeneID:104645797;Name=LOC10'
    mrna_string = 'NC_015438.2\tBestRefSeq\tmRNA\t13024\t15024\t.\t+\t.\tID=rna0;Parent=gene0;Dbxref=GeneID:'
    exon_string = 'NC_015438.2\tGnomon\texon\t4343\t4809\t.\t+\t.\tID=id1;Parent=rna0;Dbxref=GeneID:104645797'
    gene_entry = gffhelper.GFFObject(gene_string)
    handler = gffimporter.SuperLocusHandler()
    handler.gen_data_from_gffentry(gene_entry)
    handler.data.sequence_info = slice
    print(sliceh.gffid_to_coords.keys())
    print(sliceh._gff_seq_ids)
    sess.add(handler.data)
    sess.commit()
    assert handler.data.given_id == 'gene0'
    assert handler.data.type.value == 'gene'

    mrna_entry = gffhelper.GFFObject(mrna_string)
    mrna_handler = gffimporter.TranscribedHandler()
    mrna_handler.gen_data_from_gffentry(mrna_entry, super_locus=handler.data)
    piece_handler = gffimporter.TranscribedPieceHandler()
    piece_handler.gen_data_from_gffentry(mrna_entry, super_locus=handler.data, transcribed=mrna_handler.data)

    sess.add_all([mrna_handler.data, piece_handler.data])
    sess.commit()
    assert mrna_handler.data.given_id == 'rna0'
    assert mrna_handler.data.type.value == 'mRNA'
    assert mrna_handler.data.super_locus is handler.data

    exon_entry = gffhelper.GFFObject(exon_string)
    controller.clean_entry(exon_entry)
    exon_handler = gffimporter.FeatureHandler()
    exon_handler.process_gffentry(exon_entry, super_locus=handler.data, transcribed_pieces=[piece_handler.data],
                                  coordinates=coors)

    d = exon_handler.data
    s = """
    seqid {} {}
    start {} {}
    is_plus_strand {} {}
    score {} {}
    source {} {}
    phase {} {}
    given_id {} {}""".format(d.coordinates.seqid, type(d.coordinates.seqid),
                             d.position, type(d.position),
                             d.is_plus_strand, type(d.is_plus_strand),
                             d.score, type(d.score),
                             d.source, type(d.source),
                             d.phase, type(d.phase),
                             d.given_id, type(d.given_id))
    print(s)
    sess.add(exon_handler.data)
    sess.commit()

    assert exon_handler.gffentry.start == 4343
    assert exon_handler.data.is_plus_strand
    assert exon_handler.data.score is None
    assert exon_handler.data.coordinates.seqid == 'NC_015438.2'
    assert exon_handler.data.type.value == 'exon'
    assert exon_handler.data.super_locus is handler.data
    assert piece_handler.data in exon_handler.data.transcribed_pieces
    assert exon_handler.data.translateds == []
    assert exon_handler.data.transcribed_pieces[0].transcribed == mrna_handler.data


def test_data_from_cds_gffentry():
    s = "NC_015447.2\tGnomon\tCDS\t5748\t5840\t.\t-\t0\tID=cds28210;Parent=rna33721;Dbxref=GeneID:101263940,Genbank:" \
        "XP_004248424.1;Name=XP_004248424.1;gbkey=CDS;gene=LOC101263940;product=protein IQ-DOMAIN 14-like;" \
        "protein_id=XP_004248424.1"
    cds_entry = gffhelper.GFFObject(s)
    controller = gffimporter.ImportControl(database_path='sqlite:///:memory:', err_path=None)
    slic, slich = setup_data_handler(gffimporter.SequenceInfoHandler, orm.SequenceInfo)
    coords = orm.Coordinates(sequence_info=slic, seqid='dummy')
    controller.clean_entry(cds_entry)
    handler = gffimporter.FeatureHandler()
    handler.gen_data_from_gffentry(cds_entry)
    print([x.value for x in types.OnSequence])
    controller.session.add(handler.data)
    controller.session.commit()
    assert not handler.data.is_plus_strand
    assert handler.data.type.value == 'CDS'
    assert handler.data.phase == 0
    assert handler.data.score is None


def setup_testable_super_loci(db_path='sqlite:///:memory:'):
    controller = gffimporter.ImportControl(err_path='/dev/null', database_path=db_path)
    controller.mk_session()
    controller.add_sequences('testdata/dummyloci.fa')
    controller.add_gff('testdata/dummyloci.gff3', clean=False)
    return controller.super_loci[0], controller


def test_organize_and_split_features():
    sl, _ = setup_testable_super_loci()
    transcript_full = [x for x in sl.transcribed_handlers if x.data.given_id == 'y']
    print([x.data.given_id for x in sl.transcribed_handlers])
    assert len(transcript_full) == 1
    transcript_full = transcript_full[0]
    transcript_interpreter = gffimporter.TranscriptInterpreter(transcript_full)
    ordered_features = transcript_interpreter.organize_and_split_features()
    ordered_features = list(ordered_features)
    for i in [0, 4]:
        assert len(ordered_features[i]) == 1
        assert 'CDS' not in [x.data.data.type.value for x in ordered_features[i]]
    for i in [1, 2, 3]:
        assert len(ordered_features[i]) == 2
        assert 'CDS' in [x.data.data.type.value for x in ordered_features[i]]

    transcript_short = [x for x in sl.transcribed_handlers if x.data.given_id == 'z'][0]
    transcript_interpreter = gffimporter.TranscriptInterpreter(transcript_short)
    ordered_features = transcript_interpreter.organize_and_split_features()
    ordered_features = list(ordered_features)
    assert len(ordered_features) == 1
    assert len(ordered_features[0]) == 2


def test_possible_types():
    cds = types.OnSequence.CDS.name
    five_prime = types.OnSequence.five_prime_UTR.name
    three_prime = types.OnSequence.three_prime_UTR.name

    sl, _ = setup_testable_super_loci()
    transcript_full = [x for x in sl.transcribed_handlers if x.data.given_id == 'y']
    transcript_full = transcript_full[0]
    transcript_interpreter = gffimporter.TranscriptInterpreter(transcript_full)
    ordered_features = transcript_interpreter.intervals_5to3(plus_strand=True)
    ordered_features = list(ordered_features)
    pt = transcript_interpreter.possible_types(ordered_features[0])
    assert set(pt) == {five_prime, three_prime}
    pt = transcript_interpreter.possible_types(ordered_features[1])
    assert set(pt) == {cds}
    pt = transcript_interpreter.possible_types(ordered_features[-1])
    assert set(pt) == {five_prime, three_prime}


def test_import_seqinfo():
    controller = gffimporter.ImportControl(database_path='sqlite:///:memory:')
    controller.mk_session()
    seq_path = 'testdata/dummyloci.fa'
    controller.add_sequences(seq_path)
    coors = controller.sequence_info.data.coordinates
    assert len(coors) == 1
    assert coors[0].seqid == '1'
    assert coors[0].start == 0
    assert coors[0].end == 405
    assert coors[0].sha1 == 'dc6f3ba2b0c08f7d08053837b810f86cbaa06f38'  # sha1 for 'N' * 405


def test_fullcopy():
    sess = mk_session()
    sl, slh = setup_data_handler(api.SuperLocusHandler, orm.SuperLocus)
    scribed, scribedh = setup_data_handler(api.TranscribedHandler, orm.Transcribed, super_locus=sl)
    scribedpiece, scribedpieceh = setup_data_handler(api.TranscribedPieceHandler,
                                                     orm.TranscribedPiece, transcribed=scribed,
                                                     super_locus=sl, given_id='soup', )
    f, fh = setup_data_handler(api.FeatureHandler, orm.Feature, super_locus=sl,
                               transcribed_pieces=[scribedpiece], position=13)
    sess.add_all([scribedpiece, f])
    sess.commit()
    tdict = orm.Transcribed.__dict__
    print(tdict.keys())
    for key in tdict.keys():
        print('{} {}'.format(key, type(tdict[key])))

    # try copying feature
    new, newh = setup_data_handler(api.FeatureHandler, orm.Feature)
    fh.fax_all_attrs_to_another(newh)
    sess.commit()
    assert new.position == 13
    assert set(scribedpiece.features) == {f, new}
    assert new.super_locus == sl
    assert new is not f
    assert new.id != f.id

    # try copying most of transcribed things to translated
    slated, slatedh = setup_data_handler(api.TranslatedHandler, orm.Translated)
    print(scribedpiece.transcribed, 'piece.transcriebd')
    scribedpieceh.fax_all_attrs_to_another(slatedh, skip_copying=None, skip_linking=None)
    sess.commit()
    assert slated.given_id == 'soup'
    assert set(slated.features) == {f, new}
    assert f.translateds == [slated]
    assert new.transcribed_pieces == [scribedpiece]



def test_transcript_interpreter():
    sl, controller = setup_testable_super_loci()
    transcript = [x for x in sl.data.transcribeds if x.given_id == 'y'][0]
    # change so that there are implicit UTRs
    t_interp = gffimporter.TranscriptInterpreter(transcript.handler)
    t_interp.decode_raw_features()
    controller.session.commit()
    # has all standard features
    types_out = set([x.data.type.value for x in t_interp.clean_features])
    assert types_out == {types.CODING,
                         types.TRANSCRIBED,
                         types.INTRON}
    bearings_out = set([x.data.bearing.value for x in t_interp.clean_features])
    assert bearings_out == {types.START, types.END}

    assert t_interp.clean_features[-1].data.position == 400
    assert t_interp.clean_features[0].data.position == 0


def test_transcript_get_first():
    # plus strand
    sl, controller = setup_testable_super_loci()
    transcript = [x for x in sl.data.transcribeds if x.given_id == 'y'][0]
    t_interp = gffimporter.TranscriptInterpreter(transcript.handler)
    i0 = t_interp.intervals_5to3(plus_strand=True)[0]
    t_interp.interpret_first_pos(i0)
    features = t_interp.clean_features
    status = t_interp.status
    assert len(features) == 1
    f0 = features[0]
    print(f0)
    print(status.__dict__)
    print(i0[0].data.data.is_plus_strand)
    assert f0.data.position == 0
    assert status.is_5p_utr()
    assert f0.data.phase is None
    assert f0.data.is_plus_strand

    # minus strand
    sl, controller = setup_testable_super_loci()
    transcript = [x for x in sl.data.transcribeds if x.given_id == 'y'][0]
    for feature in sl.data.features:  # force minus strand
        feature.is_plus_strand = False

    # new transcript interpreter so the clean features reset
    t_interp = gffimporter.TranscriptInterpreter(transcript.handler)
    i0 = t_interp.intervals_5to3(plus_strand=False)[0]
    t_interp.interpret_first_pos(i0, plus_strand=False)
    features = t_interp.clean_features
    status = t_interp.status
    assert len(features) == 1
    f0 = features[0]
    print(f0)
    print(status)
    print(i0[0].data.data.is_plus_strand)
    print(f0.data.type)
    assert f0.data.position == 399
    assert status.is_5p_utr()
    assert f0.data.phase is None
    assert not f0.data.is_plus_strand

    # test without UTR (x doesn't have last exon, and therefore will end in CDS); remember, flipped it to minus strand
    transcript = [x for x in sl.data.transcribeds if x.given_id == 'x'][0]
    t_interp = gffimporter.TranscriptInterpreter(transcript.handler)
    i0 = t_interp.intervals_5to3(plus_strand=False)[0]
    t_interp.interpret_first_pos(i0, plus_strand=False)
    features = t_interp.clean_features
    status = t_interp.status
    assert len(features) == 4
    f_err_open = features[0]
    f_err_close = features[1]
    f_status_coding = features[2]
    f_status_transcribed = features[3]
    print(f_err_open, f_err_open.data)
    print(status)
    print(i0)
    # should get in_translated_region instead of a start codon
    assert f_status_coding.data.position == 119
    assert f_status_coding.data.type == types.CODING
    assert f_status_coding.data.bearing == types.OPEN_STATUS
    assert not f_status_coding.data.is_plus_strand
    # and should get accompanying in raw transcript
    assert f_status_transcribed.data.type == types.TRANSCRIBED
    assert f_status_coding.data.bearing == types.OPEN_STATUS
    # region beyond exon should be marked erroneous
    assert not f_err_close.data.is_plus_strand and not f_err_open.data.is_plus_strand
    assert f_err_close.data.position == 118  # so that err overlaps 1bp with the coding status checked above
    assert f_err_open.data.position == 404
    assert f_err_open.data.type == types.ERROR
    assert f_err_open.data.bearing == types.START
    assert f_err_close.data.type == types.ERROR
    assert f_err_close.data.bearing == types.END
    assert status.is_coding()
    assert status.seen_start
    assert status.genic


def test_transcript_transition_from_5p_to_end():
    sl, controller = setup_testable_super_loci()
    transcript = [x for x in sl.data.transcribeds if x.given_id == 'y'][0]
    t_interp = gffimporter.TranscriptInterpreter(transcript.handler)
    ivals_sets = t_interp.intervals_5to3(plus_strand=True)
    t_interp.interpret_first_pos(ivals_sets[0])
    # hit start codon
    t_interp.interpret_transition(ivals_before=ivals_sets[0], ivals_after=ivals_sets[1], plus_strand=True)
    features = t_interp.clean_features
    assert features[-1].data.type == types.CODING
    assert features[-1].data.bearing == types.START
    assert features[-1].data.position == 10
    # hit splice site
    t_interp.interpret_transition(ivals_before=ivals_sets[1], ivals_after=ivals_sets[2], plus_strand=True)
    assert features[-1].data.type == types.INTRON
    assert features[-1].data.bearing == types.END
    assert features[-2].data.type == types.INTRON
    assert features[-2].data.bearing == types.START
    assert features[-2].data.position == 100  # splice from
    assert features[-1].data.position == 110  # splice to
    assert t_interp.status.is_coding()
    # hit splice site
    t_interp.interpret_transition(ivals_before=ivals_sets[2], ivals_after=ivals_sets[3], plus_strand=True)
    assert features[-1].data.type == types.INTRON
    assert features[-1].data.bearing == types.END
    assert features[-2].data.type == types.INTRON
    assert features[-2].data.bearing == types.START
    assert features[-2].data.position == 120  # splice from
    assert features[-1].data.position == 200  # splice to
    # hit stop codon
    t_interp.interpret_transition(ivals_before=ivals_sets[3], ivals_after=ivals_sets[4], plus_strand=True)
    assert features[-1].data.type == types.CODING
    assert features[-1].data.bearing == types.END
    assert features[-1].data.position == 300
    # hit transcription termination site
    t_interp.interpret_last_pos(ivals_sets[4], plus_strand=True)
    assert features[-1].data.type == types.TRANSCRIBED
    assert features[-1].data.bearing == types.END
    assert features[-1].data.position == 400


def test_non_coding_transitions():
    sl, controller = setup_testable_super_loci()
    transcript = [x for x in sl.data.transcribeds if x.given_id == 'z'][0]
    piece = transcript.transcribed_pieces[0]
    t_interp = gffimporter.TranscriptInterpreter(transcript.handler)
    # get single-exon no-CDS transcript
    cds = [x for x in transcript.transcribed_pieces[0].features if x.type.value == types.CDS][0]
    piece.handler.de_link(cds.handler)
    print(transcript)
    ivals_sets = t_interp.intervals_5to3(plus_strand=True)
    assert len(ivals_sets) == 1
    t_interp.interpret_first_pos(ivals_sets[0])
    features = t_interp.clean_features
    assert features[-1].data.type == types.TRANSCRIBED
    assert features[-1].data.bearing == types.START
    assert features[-1].data.position == 110
    t_interp.interpret_last_pos(ivals_sets[0], plus_strand=True)
    assert features[-1].data.type == types.TRANSCRIBED
    assert features[-1].data.bearing == types.END
    assert features[-1].data.position == 120
    assert len(features) == 2


def test_errors_not_lost():
    sl, controller = setup_testable_super_loci()
    s = "1\tGnomon\tgene\t20\t405\t0.\t-\t0\tID=eg_missing_children"
    gene_entry = gffhelper.GFFObject(s)

    coordinates = controller.sequence_info.data.coordinates[0]
    controller.session.add(coordinates)
    sl._mark_erroneous(gene_entry, coordinates=coordinates)
    print(sl.data.transcribeds, len(sl.data.transcribeds), '...sl transcribeds')
    feature_eh, feature_e2h = sl.feature_handlers[-2:]

    print('what features did we start with::?')
    for feature in sl.data.features:
        print(feature)
        controller.session.add(feature)
    controller.session.commit()

    sl.check_and_fix_structure(sess=controller.session, coordinates=coordinates)
    for feature in sl.data.features:
        controller.session.add(feature)
    controller.session.commit()
    print('---and what features did we leave?---')
    for feature in sl.data.features:
        print(feature)
    print(feature_eh.delete_me)
    print(str(feature_eh.data), 'hello...')
    # note, you probably get sqlalchemy.orm.exc.DetachedInstanceError before failing on AssertionError below
    assert feature_eh.data in sl.data.features
    assert feature_e2h.data in sl.data.features


def test_setup_proteins():
    sl, controller = setup_testable_super_loci()
    transcript = [x for x in sl.data.transcribeds if x.given_id == 'y'][0]
    t_interp = gffimporter.TranscriptInterpreter(transcript.handler)
    print(t_interp.proteins)
    assert len(t_interp.proteins.keys()) == 1
    protein = t_interp.proteins['y.p'].data
    assert transcript in protein.transcribeds
    assert protein in transcript.translateds
    assert protein.given_id == 'y.p'
    controller.session.commit()


def test_mv_features_to_prot():
    sl, controller = setup_testable_super_loci()
    controller.session.commit()
    transcript = [x for x in sl.data.transcribeds if x.given_id == 'y'][0]
    t_interp = gffimporter.TranscriptInterpreter(transcript.handler)
    protein = t_interp.proteins['y.p'].data
    proteins = controller.session.query(orm.Translated).all()
    print(proteins, '\n^-- proteins')
    print([p.id for p in proteins])
    print(protein.id, 'protein id at start')
    conn = controller.engine.connect()
    # todo, going from TUESDAY, whyyyyyyyyy?
    conn.execute(orm.association_translateds_to_features.insert(), [{'translated_id': 1, 'feature_id': 14}])
    t_interp.decode_raw_features()
    t_interp.mv_coding_features_to_proteins(controller.feature2protein_to_add)
    controller.session.commit()
    # grab protein again jic
    protein = controller.session.query(orm.Translated).filter(orm.Translated.given_id == 'y.p').all()
    assert len(protein) == 1
    protein = protein[0]
    print(protein.id)
    print(controller.session.query(orm.association_translateds_to_features).all())
    assert len(protein.features) == 2  # start and stop codon
    assert set([x.type.value for x in protein.features]) == {types.CODING}
    assert set([x.bearing.value for x in protein.features]) == {types.START, types.END}


def test_check_and_fix_structure():
    rel_path = 'testdata/dummyloci_annotations.sqlitedb'  # so we save a copy of the cleaned up loci once per test run
    db_path = 'sqlite:///{}'.format(rel_path)
    if os.path.exists(rel_path):
        os.remove(rel_path)
    sl, controller = setup_testable_super_loci(db_path)
    coordinates = controller.sequence_info.data.coordinates[0]

    sl.check_and_fix_structure(controller.session, coordinates=coordinates, controller=controller)
    # check handling of nice transcript
    transcript = [x for x in sl.data.transcribeds if x.given_id == 'y'][0]
    #protein = [x for x in sl.data.translateds if x.given_id == 'y.p'][0]
    protein = controller.session.query(orm.Translated).filter(orm.Translated.given_id == 'y.p').all()

    assert len(protein) == 1
    protein = protein[0]
    print(protein.id)
    print(controller.session.query(orm.association_translateds_to_features).all())
    # check we get a protein with start and stop codon for the nice transcript
    assert len(protein.features) == 2  # start and stop codon
    assert set([x.type.value for x in protein.features]) == {types.CODING}
    assert set([x.position for x in protein.features]) == {10, 300}
    assert set([x.bearing.value for x in protein.features]) == {types.START, types.END}
    # check we get a transcript with tss, 2x(dss, ass), and tts (+ start & stop codons)
    piece = transcript.handler.one_piece().data
    assert len(piece.features) == 8
    assert set([x.type.value for x in piece.features]) == {types.TRANSCRIBED,
                                                           types.INTRON,
                                                           types.CODING,
                                                           }
    assert set([x.bearing.value for x in piece.features]) == {types.START, types.END}
    # check handling of truncated transcript
    transcript = [x for x in sl.data.transcribeds if x.given_id == 'x'][0]
    piece = transcript.handler.one_piece().data
    protein = [x for x in sl.data.translateds if x.given_id == 'x.p'][0]
    print(protein.features)
    assert len(protein.features) == 2
    assert set([x.type.value for x in protein.features]) == {types.CODING}
    assert set([x.position for x in protein.features]) == {10, 120}
    assert set([x.bearing.value for x in protein.features]) == {types.START, types.CLOSE_STATUS}

    assert len(piece.features) == 8
    assert set([x.type.value for x in piece.features]) == {types.TRANSCRIBED, types.INTRON,
                                                           types.ERROR, types.CODING}
    coding_fs = [x for x in piece.features if x.type.value == types.CODING]
    assert len(coding_fs) == 2
    assert set([x.bearing.value for x in coding_fs]) == {types.START, types.CLOSE_STATUS}

    transcribed_fs = [x for x in piece.features if x.type.value == types.TRANSCRIBED]
    assert len(transcribed_fs) == 2
    assert set([x.bearing.value for x in transcribed_fs]) == {types.START, types.CLOSE_STATUS}
    assert set([x.position for x in transcribed_fs]) == {0, 120}

    assert len(sl.data.translateds) == 3
    controller.session.commit()


def test_erroneous_splice():
    db_path = 'sqlite:///:memory:'

    sl, controller = setup_testable_super_loci(db_path)
    sess = controller.session
    # get target transcript
    transcript = [x for x in sl.data.transcribeds if x.given_id == 'x'][0]
    # fish out "first exon" features and extend so intron is of -length
    f0 = sess.query(orm.Feature).filter(orm.Feature.given_id == 'ftr000000').first()
    f1 = sess.query(orm.Feature).filter(orm.Feature.given_id == 'ftr000001').first()
    f0.handler.gffentry.end = f1.handler.gffentry.end = 115

    ti = gffimporter.TranscriptInterpreter(transcript.handler)
    ti.decode_raw_features()
    clean_datas = [x.data for x in ti.clean_features]
    # TSS, start codon, 2x error splice, 2x error splice, 2x error no stop
    print('---\n'.join([str(x) for x in clean_datas]))
    assert len(clean_datas) == 10

    assert len([x for x in clean_datas if x.type == types.ERROR]) == 6
    # make sure splice error covers whole exon-intron-exon region
    assert clean_datas[2].type == types.ERROR
    assert clean_datas[2].bearing == types.START
    assert clean_datas[2].position == 10
    assert clean_datas[3].type == types.ERROR
    assert clean_datas[3].bearing == types.END
    assert clean_datas[3].position == 120
    sess.commit()


def test_gff_gen():
    controller = gffimporter.ImportControl(database_path='sqlite:///:memory:')
    x = list(controller.gff_gen('testdata/testerSl.gff3'))
    assert len(x) == 103
    assert x[0].type == 'region'
    assert x[-1].type == 'CDS'


def test_gff_useful_gen():
    controller = gffimporter.ImportControl(database_path='sqlite:///:memory:')
    x = list(controller.useful_gff_entries('testdata/testerSl.gff3'))
    assert len(x) == 100  # should drop the region entry
    assert x[0].type == 'gene'
    assert x[-1].type == 'CDS'


def test_gff_grouper():
    controller = gffimporter.ImportControl(database_path='sqlite:///:memory:')
    x = list(controller.group_gff_by_gene('testdata/testerSl.gff3'))
    assert len(x) == 5
    for group in x:
        assert group[0].type == 'gene'

# section: types
def test_enum_non_inheritance():
    allknown = [x.name for x in list(types.AllKnown)]
    allnice = [x.name for x in list(types.AllKeepable)]
    # check that some random bits made it in to all
    assert 'error' in allknown
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

