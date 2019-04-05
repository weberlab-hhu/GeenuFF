from . import types

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Table, Column, Integer, ForeignKey, String, Enum, CheckConstraint, UniqueConstraint, Boolean, Float
from sqlalchemy.orm import relationship

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


class Coordinate(Base):
    __tablename__ = 'coordinate'

    id = Column(Integer, primary_key=True)
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
        {})

    def __repr__(self):
        return '<Coordinate {}, {}:{}-{}>'.format(self.id, self.seqid, self.start, self.end)


class SuperLocusAliases(Base):
    __tablename__ = 'super_locus_aliases'

    id = Column(Integer, primary_key=True)
    alias = Column(String)
    super_locus_id = Column(Integer, ForeignKey('super_locus.id'), nullable=False)
    super_locus = relationship('SuperLocus', back_populates='aliases')


class SuperLocus(Base):
    __tablename__ = 'super_locus'
    # normally a loci, some times a short list of loci for "trans splicing"
    # this will define a group of exons that can possibly be made into transcripts
    # AKA this if you have to go searching through a graph for parents/children, at least said graph will have
    # a max size defined at SuperLoci

    id = Column(Integer, primary_key=True)
    given_id = Column(String)
    type = Column(Enum(types.SuperLocusAll))
    # things SuperLocus can have a lot of
    aliases = relationship('SuperLocusAliases', back_populates='super_locus')
    transcribeds = relationship('Transcribed', back_populates='super_locus')
    translateds = relationship('Translated', back_populates='super_locus')


association_transcribed_piece_to_feature = Table('association_transcribed_piece_to_feature', Base.metadata,  # todo, rename
    Column('transcribed_piece_id', Integer, ForeignKey('transcribed_piece.id'), nullable=False),
    Column('feature_id', Integer, ForeignKey('feature.id'), nullable=False)
)


association_translated_to_feature = Table('association_translated_to_feature', Base.metadata,
    Column('translated_id', Integer, ForeignKey('translated.id'), nullable=False),
    Column('feature_id', Integer, ForeignKey('feature.id'), nullable=False)
)


class Transcribed(Base):
    __tablename__ = 'transcribed'

    id = Column(Integer, primary_key=True)
    given_id = Column(String)

    type = Column(Enum(types.TranscriptLevelAll))

    super_locus_id = Column(Integer, ForeignKey('super_locus.id'), nullable=False)
    super_locus = relationship('SuperLocus', back_populates='transcribeds')

    transcribed_pieces = relationship('TranscribedPiece', back_populates='transcribed')

    def __repr__(self):
        return '<Transcribed, {}, "{}" of type {}, with {} pieces>'.format(self.id, self.given_id, self.type,
                                                                           len(self.transcribed_pieces))


class TranscribedPiece(Base):
    __tablename__ = 'transcribed_piece'

    id = Column(Integer, primary_key=True)
    given_id = Column(String)

    position = Column(Integer, nullable=False)

    transcribed_id = Column(Integer, ForeignKey('transcribed.id'), nullable=False)
    transcribed = relationship('Transcribed', back_populates='transcribed_pieces')

    features = relationship('Feature', secondary=association_transcribed_piece_to_feature,
                            back_populates='transcribed_pieces')

    __table_args__ = (UniqueConstraint('transcribed_id', 'position'),)

    def __repr__(self):
        return "<TranscribedPiece, {}: with features {}>".format(
            self.id, [(x.id, x.position, x.given_id) for x in self.features])


class Translated(Base):
    __tablename__ = 'translated'

    id = Column(Integer, primary_key=True)
    given_id = Column(String)
    # type can only be 'protein' so far as I know..., so skipping
    super_locus_id = Column(Integer, ForeignKey('super_locus.id'), nullable=False)
    super_locus = relationship('SuperLocus', back_populates='translateds')

    features = relationship('Feature', secondary=association_translated_to_feature,
                            back_populates='translateds')



class Feature(Base):
    __tablename__ = 'feature'
    # basic attributes
    id = Column(Integer, primary_key=True)
    given_id = Column(String)
    type = Column(Enum(types.OnSequence))
    bearing = Column(Enum(types.Bearings))
    position = Column(Integer)
    is_plus_strand = Column(Boolean)
    score = Column(Float)
    source = Column(String)
    phase = Column(Integer)

    # for differentiating from subclass entries
    subtype = Column(String(20))

    #seqid = Column(String)
    # any piece of coordinate always has just one seqid
    coordinate_id = Column(Integer, ForeignKey('coordinate.id'), nullable=False)
    coordinate = relationship('Coordinate', back_populates='features')

    # relations
    transcribed_pieces = relationship('TranscribedPiece',
                                      secondary=association_transcribed_piece_to_feature,
                                      back_populates='features')

    translateds = relationship('Translated',
                               secondary=association_translated_to_feature,
                               back_populates='features')

    __table_args__ = (
        CheckConstraint(position >= -1, name='check_start_1plus'),
        CheckConstraint(phase >= 0, name='check_phase_not_negative'),
        CheckConstraint(phase < 3, name='check_phase_less_three'),
    )

    __mapper_args__ = {
        'polymorphic_on': subtype,
        'polymorphic_identity': 'general'
    }

    def __repr__(self):
        s = '<{py_type}, {pk}: {givenid} of type: {type} ({bearing}) @{position} on {coor}, is_plus: {plus}, ' \
            'phase: {phase}>'.format(
                pk=self.id, bearing=self.bearing,
                type=self.type, position=self.position, coor=self.coordinate, plus=self.is_plus_strand,
                phase=self.phase, givenid=self.given_id, py_type=type(self)
            )
        return s

    def cmp_key(self):  # todo, pos_cmp & full_cmp
        return self.coordinate.seqid, self.is_plus_strand, self.position, self.type

    def pos_cmp_key(self):
        return self.coordinate.seqid, self.is_plus_strand, self.position


