from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Table, Column, Integer, ForeignKey, String, Enum, CheckConstraint, UniqueConstraint, Boolean, Float
from sqlalchemy.orm import relationship

from . import types

# setup classes for data holding
Base = declarative_base()


class Genome(Base):
    __tablename__ = 'genome'

    id = Column(Integer, primary_key=True)
    species = Column(String)
    accession = Column(String)
    version = Column(String)
    acquired_from = Column(String)
    coordinates = relationship("Coordinate", back_populates="genome")

    __table_args__ = (
        UniqueConstraint('species', 'version', name='unique_genome_species_and_version'),
    )

    def __repr__(self):
        return '<Genome {}, species: {}>'.format(self.id, self.species)


class Coordinate(Base):
    __tablename__ = 'coordinate'

    id = Column(Integer, primary_key=True)
    sequence = Column(String)
    length = Column(Integer)
    seqid = Column(String, nullable=False)
    sha1 = Column(String)
    genome_id = Column(Integer, ForeignKey('genome.id'), nullable=False, index=True)
    genome = relationship('Genome', back_populates='coordinates')

    features = relationship('Feature', back_populates='coordinate')

    __table_args__ = (
        UniqueConstraint('genome_id', 'seqid', name='unique_coords_per_genome'),
        CheckConstraint('length >= 0', name='check_positive_length'),
    )

    def __repr__(self):
        return '<Coordinate {}, seqid: {}, len: {}, genome_id: {}, genome: {}>'.format(self.id, self.seqid, self.length,
                                                                                       self.genome_id,
                                                                                       self.genome.species)


class SuperLocus(Base):
    __tablename__ = 'super_locus'
    # normally a loci, some times a short list of loci for "trans splicing"
    # this will define a group of exons that can possibly be made into transcripts
    # AKA this if you have to go searching through a graph for parents/children, at least said graph will have
    # a max size defined at SuperLoci

    id = Column(Integer, primary_key=True)
    given_name = Column(String)
    aliases = Column(String)
    type = Column(Enum(types.SuperLocusAll), index=True)
    # things SuperLocus can have a lot of
    transcripts = relationship('Transcript', back_populates='super_locus')
    proteins = relationship('Protein', back_populates='super_locus')

    def __repr__(self):
        return '<SuperLocus {}, given_name: \'{}\', type: {}>'.format(self.id, self.given_name,
                                                                      self.type.value)


association_transcript_piece_to_feature = Table('association_transcript_piece_to_feature', Base.metadata,
    Column('transcript_piece_id', Integer, ForeignKey('transcript_piece.id'), nullable=False, index=True),
    Column('feature_id', Integer, ForeignKey('feature.id'), nullable=False, index=True)
)


association_protein_to_feature = Table('association_protein_to_feature', Base.metadata,
    Column('protein_id', Integer, ForeignKey('protein.id'), nullable=False, index=True),
    Column('feature_id', Integer, ForeignKey('feature.id'), nullable=False, index=True)
)

association_transcript_to_protein = Table('association_transcript_to_protein', Base.metadata,
    Column('transcript_id', Integer, ForeignKey('transcript.id'), nullable=False, index=True),
    Column('protein_id', Integer, ForeignKey('protein.id'), nullable=False, index=True)
)


class Transcript(Base):
    __tablename__ = 'transcript'

    id = Column(Integer, primary_key=True)
    given_name = Column(String)

    type = Column(Enum(types.TranscriptLevel), index=True)
    longest = Column(Boolean, index=True)

    super_locus_id = Column(Integer, ForeignKey('super_locus.id'), nullable=False, index=True)
    super_locus = relationship('SuperLocus', back_populates='transcripts')

    transcript_pieces = relationship('TranscriptPiece', back_populates='transcript')

    proteins = relationship('Protein', secondary=association_transcript_to_protein,
                            back_populates='transcripts')

    def __repr__(self):
        return '<Transcript, {}, "{}" of type {}, with {} pieces>'.format(self.id,
                                                                          self.given_name,
                                                                          self.type,
                                                                          len(self.transcript_pieces))


class TranscriptPiece(Base):
    __tablename__ = 'transcript_piece'

    id = Column(Integer, primary_key=True)
    given_name = Column(String)

    position = Column(Integer, nullable=False)

    transcript_id = Column(Integer, ForeignKey('transcript.id'), nullable=False, index=True)
    transcript = relationship('Transcript', back_populates='transcript_pieces')

    features = relationship('Feature', secondary=association_transcript_piece_to_feature,
                            back_populates='transcript_pieces')

    __table_args__ = (
        UniqueConstraint('transcript_id', 'position', name='unique_positions_per_piece'),
    )

    def __repr__(self):
        return ('<TranscriptPiece, {}: for transcript {} '
                'in position {}>').format(self.id, self.transcript_id, self.position)


class Protein(Base):
    __tablename__ = 'protein'

    id = Column(Integer, primary_key=True)
    given_name = Column(String)
    # type can only be 'protein' so far as I know..., so skipping
    super_locus_id = Column(Integer, ForeignKey('super_locus.id'), nullable=False)
    super_locus = relationship('SuperLocus', back_populates='proteins')

    features = relationship('Feature', secondary=association_protein_to_feature,
                            back_populates='proteins')

    transcripts = relationship('Transcript', secondary=association_transcript_to_protein,
                               back_populates='proteins')

    def __repr__(self):
        return '<Protein {}, given_name: \'{}\', super_locus_id: {}>'.format(self.id,
                                                                             self.given_name,
                                                                             self.super_locus_id)


class Feature(Base):
    __tablename__ = 'feature'
    # basic attributes
    id = Column(Integer, primary_key=True)
    given_name = Column(String)
    type = Column(Enum(types.GeenuffFeature), index=True)

    start = Column(Integer, nullable=False)
    start_is_biological_start = Column(Boolean, nullable=False)
    end = Column(Integer, nullable=False)
    end_is_biological_end = Column(Boolean, nullable=False)

    is_plus_strand = Column(Boolean, nullable=False)
    score = Column(Float)
    source = Column(String)
    phase = Column(Integer)

    # any piece of coordinate always has just one seqid
    coordinate_id = Column(Integer, ForeignKey('coordinate.id'), nullable=False, index=True)
    coordinate = relationship('Coordinate', back_populates='features')

    # relations
    transcript_pieces = relationship('TranscriptPiece',
                                      secondary=association_transcript_piece_to_feature,
                                      back_populates='features')

    proteins = relationship('Protein',
                               secondary=association_protein_to_feature,
                               back_populates='features')

    __table_args__ = (
        UniqueConstraint('coordinate_id', 'type', 'start', 'end', 'is_plus_strand', 'given_name',
                         name='unique_feature'),
        CheckConstraint('start >= 0 and end >= -1 and phase >= 0 and phase < 3',
                        name='check_start_end_phase'),
        # if start_is_biological_start is True, phase has to be 0
        CheckConstraint('phase is NULL or (not start_is_biological_start or phase = 0)',
                        name='check_phase_bio_start'),
        # check start/end order depending on is_plus_strand
        CheckConstraint('(is_plus_strand and start <= end) or (not is_plus_strand and end <= start)',
                        name='start_end_order')
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
