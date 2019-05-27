from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Table, Column, Integer, ForeignKey, String, Enum, CheckConstraint, UniqueConstraint, Boolean, Float
from sqlalchemy.orm import relationship

from . import types

# setup classes for data holding
Base = declarative_base()


class Genome(Base):
    __tablename__ = 'genome'

    # data
    id = Column(Integer, primary_key=True)
    species = Column(String)
    accession = Column(String)
    version = Column(String)
    acquired_from = Column(String)
    coordinates = relationship("Coordinate", back_populates="genome")

    def __repr__(self):
        return '<Genome {}, species: {}>'.format(self.id, self.species)

class Coordinate(Base):
    __tablename__ = 'coordinate'

    id = Column(Integer, primary_key=True)
    sequence = Column(String)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    seqid = Column(String, nullable=False)
    sha1 = Column(String)
    genome_id = Column(Integer, ForeignKey('genome.id'), nullable=False)
    genome = relationship('Genome', back_populates='coordinates')

    features = relationship('Feature', back_populates='coordinate')

    __table_args__ = (
        CheckConstraint(start >= 0, name='check_start_1plus'),
        CheckConstraint(end > start, name='check_end_gr_start'),
    )

    def __repr__(self):
        return '<Coordinate {}, {}:{}--{}>'.format(self.id, self.seqid, self.start, self.end)


class SuperLocus(Base):
    __tablename__ = 'super_locus'
    # normally a loci, some times a short list of loci for "trans splicing"
    # this will define a group of exons that can possibly be made into transcripts
    # AKA this if you have to go searching through a graph for parents/children, at least said graph will have
    # a max size defined at SuperLoci

    id = Column(Integer, primary_key=True)
    given_name = Column(String)
    aliases = Column(String)
    type = Column(Enum(types.SuperLocusAll))
    # things SuperLocus can have a lot of
    transcribeds = relationship('Transcript', back_populates='super_locus')
    translateds = relationship('Protein', back_populates='super_locus')

    def __repr__(self):
        return '<SuperLocus {}, given_name: \'{}\', type: {}>'.format(self.id, self.given_name,
                                                                  self.type.value)


association_transcribed_piece_to_feature = Table('association_transcribed_piece_to_feature', Base.metadata,
    Column('transcribed_piece_id', Integer, ForeignKey('transcribed_piece.id'), nullable=False),
    Column('feature_id', Integer, ForeignKey('feature.id'), nullable=False)
)


association_translated_to_feature = Table('association_translated_to_feature', Base.metadata,
    Column('translated_id', Integer, ForeignKey('translated.id'), nullable=False),
    Column('feature_id', Integer, ForeignKey('feature.id'), nullable=False)
)


class Transcript(Base):
    __tablename__ = 'transcribed'

    id = Column(Integer, primary_key=True)
    given_name = Column(String)

    type = Column(Enum(types.TranscriptLevelAll))

    super_locus_id = Column(Integer, ForeignKey('super_locus.id'), nullable=False)
    super_locus = relationship('SuperLocus', back_populates='transcribeds')

    transcribed_pieces = relationship('TranscriptPiece', back_populates='transcribed')

    def __repr__(self):
        return '<Transcript, {}, "{}" of type {}, with {} pieces>'.format(self.id, self.given_name, self.type,
                                                                           len(self.transcribed_pieces))


class TranscriptPiece(Base):
    __tablename__ = 'transcribed_piece'

    id = Column(Integer, primary_key=True)
    given_name = Column(String)

    position = Column(Integer, nullable=False)

    transcribed_id = Column(Integer, ForeignKey('transcribed.id'), nullable=False)
    transcribed = relationship('Transcript', back_populates='transcribed_pieces')

    features = relationship('Feature', secondary=association_transcribed_piece_to_feature,
                            back_populates='transcribed_pieces')

    __table_args__ = (
        UniqueConstraint('transcribed_id', 'position', name='unique_positions_per_piece'),
    )

    def __repr__(self):
        return ('<TranscriptPiece, {}: for transcribed {} '
                'in position {}>').format(self.id, self.transcribed_id, self.position)


class Protein(Base):
    __tablename__ = 'translated'

    id = Column(Integer, primary_key=True)
    given_name = Column(String)
    # type can only be 'protein' so far as I know..., so skipping
    super_locus_id = Column(Integer, ForeignKey('super_locus.id'), nullable=False)
    super_locus = relationship('SuperLocus', back_populates='translateds')

    features = relationship('Feature', secondary=association_translated_to_feature,
                            back_populates='translateds')

    def __repr__(self):
        return '<Protein {}, given_name: \'{}\', super_locus_id: {}>'.format(self.id,
                                                                                self.given_name,
                                                                                self.super_locus_id)


class Feature(Base):
    __tablename__ = 'feature'
    # basic attributes
    id = Column(Integer, primary_key=True)
    given_name = Column(String)
    type = Column(Enum(types.OnSequence))

    start = Column(Integer)
    start_is_biological_start = Column(Boolean)
    end = Column(Integer)
    end_is_biological_end = Column(Boolean)

    is_plus_strand = Column(Boolean)
    score = Column(Float)
    source = Column(String)
    phase = Column(Integer)

    # any piece of coordinate always has just one seqid
    coordinate_id = Column(Integer, ForeignKey('coordinate.id'), nullable=False)
    coordinate = relationship('Coordinate', back_populates='features')

    # relations
    transcribed_pieces = relationship('TranscriptPiece',
                                      secondary=association_transcribed_piece_to_feature,
                                      back_populates='features')

    translateds = relationship('Protein',
                               secondary=association_translated_to_feature,
                               back_populates='features')

    __table_args__ = (
        CheckConstraint(end >= -1, name='check_end_minus1plus'),
        CheckConstraint(start >= 0, name="check_start_0plus"),
        CheckConstraint(phase >= 0, name='check_phase_not_negative'),
        CheckConstraint(phase < 3, name='check_phase_less_three'),
    )

    def __repr__(self):
        start_caveat = end_caveat = ''
        if self.start_is_biological_start is None:
            start_caveat = '?'
        elif self.start_is_biological_start is False:
            start_caveat = '*'
        if self.end_is_biological_end is None:
            end_caveat = '?'
        elif self.end_is_biological_end is False:
            end_caveat = '*'
        s = ('<Feature, id: {pk}, given_name: \'{givenid}\', type: {type}, '
             '{start}{start_caveat}--{end}{end_caveat}, on {coor}, '
             'is_plus: {plus}, phase: {phase}>').format(
                pk=self.id, start_caveat=start_caveat, end_caveat=end_caveat,
                type=self.type.value, start=self.start, end=self.end, coor=self.coordinate,
                plus=self.is_plus_strand, phase=self.phase, givenid=self.given_name
            )
        return s
