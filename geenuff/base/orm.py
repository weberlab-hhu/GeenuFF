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
        UniqueConstraint('genome_id', 'seqid', name='unique_coords_per_genome'),
        CheckConstraint(start >= 0, name='check_start_1plus'),
        CheckConstraint(end > start, name='check_end_gr_start'),
    )

    def __repr__(self):
        return '<Coordinate {}, {}:{}-{}>'.format(self.id, self.seqid, self.start, self.end)


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
    transcribeds = relationship('Transcribed', back_populates='super_locus')
    translateds = relationship('Translated', back_populates='super_locus')


association_transcribed_piece_to_feature = Table('association_transcribed_piece_to_feature', Base.metadata,
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
    given_name = Column(String)

    type = Column(Enum(types.TranscriptLevelAll))

    super_locus_id = Column(Integer, ForeignKey('super_locus.id'), nullable=False)
    super_locus = relationship('SuperLocus', back_populates='transcribeds')

    transcribed_pieces = relationship('TranscribedPiece', back_populates='transcribed')

    def __repr__(self):
        return '<Transcribed, {}, "{}" of type {}, with {} pieces>'.format(self.id, self.given_name, self.type,
                                                                           len(self.transcribed_pieces))


class TranscribedPiece(Base):
    __tablename__ = 'transcribed_piece'

    id = Column(Integer, primary_key=True)
    given_name = Column(String)

    position = Column(Integer, nullable=False)

    transcribed_id = Column(Integer, ForeignKey('transcribed.id'), nullable=False)
    transcribed = relationship('Transcribed', back_populates='transcribed_pieces')

    features = relationship('Feature', secondary=association_transcribed_piece_to_feature,
                            back_populates='transcribed_pieces')

    __table_args__ = (
        UniqueConstraint('transcribed_id', 'position', name='unique_positions_per_piece'),
    )

    def __repr__(self):
        features = [(x.id, x.start, x.end, x.given_name) for x in self.features]
        return ('<TranscribedPiece, {}: for transcribed {} '
                'in position {} with features {}>').format(self.id,
                                                           self.transcribed_id,
                                                           self.position,
                                                           features)



class Translated(Base):
    __tablename__ = 'translated'

    id = Column(Integer, primary_key=True)
    given_name = Column(String)
    # type can only be 'protein' so far as I know..., so skipping
    super_locus_id = Column(Integer, ForeignKey('super_locus.id'), nullable=False)
    super_locus = relationship('SuperLocus', back_populates='translateds')

    features = relationship('Feature', secondary=association_translated_to_feature,
                            back_populates='translateds')


class Feature(Base):
    __tablename__ = 'feature'
    # basic attributes
    id = Column(Integer, primary_key=True)
    given_name = Column(String)
    type = Column(Enum(types.OnSequence))

    start = Column(Integer, nullable=False)
    start_is_biological_start = Column(Boolean, nullable=False)
    end = Column(Integer, nullable=False)
    end_is_biological_end = Column(Boolean, nullable=False)

    is_plus_strand = Column(Boolean, nullable=False)
    score = Column(Float)
    source = Column(String)
    phase = Column(Integer)

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
        UniqueConstraint('coordinate_id', 'type', 'start', 'end', 'is_plus_strand',
                         name='unique_feature'),
        CheckConstraint(end >= -1, name='check_end_minus1plus'),
        CheckConstraint(start >= 0, name="check_start_0plus"),
        CheckConstraint(phase >= 0, name='check_phase_not_negative'),
        CheckConstraint(phase < 3, name='check_phase_less_three'),
        # if start_is_biological_start is True, phase has to be 0
        CheckConstraint('not start_is_biological_start or phase = 0', name='check_phase_bio_start')
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
        s = '<{py_type}, {pk}: {givenid} of type: {type} @({start}{start_caveat} -> {end}{end_caveat}) on {coor}, ' \
            'is_plus: {plus}, phase: {phase}>'.format(
                pk=self.id, start_caveat=start_caveat, end_caveat=end_caveat,
                type=self.type, start=self.start, end=self.end, coor=self.coordinate, plus=self.is_plus_strand,
                phase=self.phase, givenid=self.given_name, py_type=type(self)
            )
        return s

    def cmp_key(self):
        pos_cmp = list(self.pos_cmp_key())
        pos_cmp.append(self.type)
        return tuple(pos_cmp)

    def pos_cmp_key(self):
        sortable_start = self.start
        sortable_end = self.end
        if not self.is_plus_strand:
            sortable_start = sortable_start * -1
            sortable_end = sortable_end * -1
        return self.coordinate.seqid, self.is_plus_strand, sortable_start, sortable_end
