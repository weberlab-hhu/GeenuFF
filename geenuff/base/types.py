import enum


# Enum makers
def join_to_enum(name, *args):
    """joins enums from args into returned enum"""
    enum_bits = []
    for cls in args:
        enum_bits += list(cls)
    out = enum.Enum(name, [(x.name, x.value) for x in enum_bits])
    return out


def make_enum(name, *args):
    """makes enum from list of strings"""
    return enum.Enum(name, [(x, x) for x in args])


# code words as variables
# SuperLocus
CODING_GENE = 'coding_gene'
NON_CODING_GENE = 'non_coding_gene'
PSEUDOGENE = 'pseudogene'
OPERON = 'operon'

# SuperLocusHistorical
GENE = 'gene'

# TranscriptLevelNice
MRNA = 'mRNA'
TRNA = 'tRNA'
RRNA = 'rRNA'
MIRNA = 'miRNA'
SNORNA = 'snoRNA'
SNRNA = 'snRNA'
SRP_RNA = 'SRP_RNA'
LNC_RNA = 'lnc_RNA'
PRE_MIRNA = 'pre_miRNA'
RNASE_MRP_RNA = 'RNase_MRP_RNA'

# TranscriptLevelInput
TRANSCRIPT = 'transcript'
PRIMARY_TRANSCRIPT = 'primary_transcript'
PSEUDOGENIC_TRANSCRIPT = 'pseudogenic_transcript'

# make all the enums
SuperLocus = make_enum('SuperLocus', CODING_GENE, NON_CODING_GENE, PSEUDOGENE, OPERON)
SuperLocusHistorical = make_enum('SuperLocusHistorical', GENE)
SuperLocusAll = join_to_enum('SuperLocusAll', SuperLocus, SuperLocusHistorical)


TranscriptLevelNice = make_enum('TranscriptLevelNice', MRNA, TRNA, RRNA,MIRNA, SNORNA, SNRNA, SRP_RNA,
                                LNC_RNA, PRE_MIRNA, RNASE_MRP_RNA)
TranscriptLevelInput = make_enum('TranscriptLevelInput', TRANSCRIPT, PRIMARY_TRANSCRIPT, PSEUDOGENIC_TRANSCRIPT)
TranscriptLevelAll = join_to_enum('TranscriptLevelAll', TranscriptLevelNice, TranscriptLevelInput)


# FEATURES
# TranscribedGeneral
#TRANSCRIPTION_START_SITE = 'TSS'  # transcription_start_site
#TRANSCRIPTION_TERMINATION_SITE = 'TTS'  # transcription_termination_site
#DONOR_SPLICE_SITE = 'DSS'  # donor_splice_site
#ACCEPTOR_SPLICE_SITE = 'ASS'  # acceptor_splice_site
# TranscribedInput
EXON = 'exon'
#INTRON = 'intron'
FIVE_PRIME_UTR = 'five_prime_UTR'
THREE_PRIME_UTR = 'three_prime_UTR'
## TranscribedStatus
#IN_RAW_TRANSCRIPT = 'in_raw_transcript'
#IN_INTRON = 'in_intron'
#IN_TRANS_INTRON = 'in_trans_intron'
## TranscribedTransSplice
#DONOR_TRANS_SPLICE_SITE = 'donor_trans_splice_site'
#ACCEPTOR_TRANS_SPLICE_SITE = 'acceptor_trans_splice_site'


# translation related features
# TranslatedInput
CDS = 'CDS'
## TranslatedGeneral
#START_CODON = 'start_codon'
#STOP_CODON = 'stop_codon'
# TranslatedStatus
#IN_TRANSLATED_REGION = 'in_translated_region'  # so between start and stop codons, could still be in an intron

# Enums

## Anything else
## errror features
#IN_ERROR = 'in_error'
#ERROR_OPEN = 'error_open'
#ERROR_CLOSE = 'error_close'
#
#ErrorFeature = make_enum('ErrorFeature', IN_ERROR, ERROR_OPEN, ERROR_CLOSE)

REGION = 'region'
CHROMOSOME = 'chromosome'
SUPERCONTIG = 'supercontig'
MATCH = 'match'
CDNA_MATCH = 'cDNA_match'

# final feature types
TRANSCRIPT_FEATURE = 'transcript_feature'
CODING = 'coding'
INTRON = 'intron'  # defined above
TRANS_INTRON = 'trans_intron'
FinalFeatures = make_enum('FinalFeatures', TRANSCRIPT_FEATURE, CODING, INTRON, TRANS_INTRON)

# error types
MISSING_UTR_5P = 'missing_utr_5p'
MISSING_UTR_3P = 'missing_utr_3p'
EMPTY_SUPER_LOCUS = 'empty_super_locus'
MISSING_START_CODON = 'missing_start_codon'
MISSING_STOP_CODON = 'missing_stop_codon'
WRONG_PHASE_5P = 'wrong_starting_phase'
MISMATCHED_PHASE_3P = 'mismatched_ending_phase'
OVERLAPPING_EXONS = 'overlapping_exons'
TOO_SHORT_INTRON = 'too_short_intron'
SL_OVERLAP_ERROR = 'super_loci_overlap_error'
Errors = make_enum('Errors', MISSING_UTR_5P, MISSING_UTR_3P, EMPTY_SUPER_LOCUS, MISSING_START_CODON,
                   MISSING_STOP_CODON, WRONG_PHASE_5P, MISMATCHED_PHASE_3P, OVERLAPPING_EXONS,
                   TOO_SHORT_INTRON, SL_OVERLAP_ERROR)

# bearings
START = 'start'
END = 'end'
OPEN_STATUS = 'open_status'
CLOSE_STATUS = 'close_status'
POINT = 'point'

Bearings = make_enum('Bearings', START, END, OPEN_STATUS, CLOSE_STATUS, POINT)

# ignorable known
IgnorableFeatures = make_enum('IgnorableFeatures', REGION, CHROMOSOME, SUPERCONTIG, MATCH, CDNA_MATCH)

# On sequence final features
KeepOnSequence = FinalFeatures
#KeepOnSequence = join_to_enum('KeepOnSequence', TranscribedNice, TranslatedNice, ErrorFeature)

# transcription related features
TranscriptGeneral = make_enum('TranscriptGeneral', TRANSCRIPT_FEATURE)
TranscriptInput = make_enum('TranscriptInput', EXON, FIVE_PRIME_UTR, THREE_PRIME_UTR)
#TranscribedStatus = make_enum('TranscribedStatus', IN_RAW_TRANSCRIPT, IN_INTRON, IN_TRANS_INTRON)
#TranscribedTransSplice = make_enum('TranscribedTransSplice', DONOR_TRANS_SPLICE_SITE, ACCEPTOR_TRANS_SPLICE_SITE)
# joining transcription related
TranscriptAll = join_to_enum('TranscriptFeature', TranscriptGeneral, TranscriptInput)
#TranscribedNice = join_to_enum('TranscribedNice', TranscribedGeneral, TranscribedStatus, TranscribedTransSplice)

ProteinInput = make_enum('ProteinInput', CDS)
ProteinGeneral = make_enum('ProteinGeneral', CODING)
#TranslatedStatus = make_enum('TranslatedStatus', IN_TRANSLATED_REGION)

ProteinAll = join_to_enum('ProteinFeatureType', ProteinInput, ProteinGeneral)
#TranslatedNice = join_to_enum('TranslatedNice', TranslatedStatus, TranslatedGeneral)


# All known features (else error on import)
OnSequence = join_to_enum('OnSequence', TranscriptInput, ProteinInput, FinalFeatures, Errors)
AllKnown = join_to_enum('AllKnown', SuperLocusAll, TranscriptLevelAll, OnSequence,
                        IgnorableFeatures)
AllKeepable = join_to_enum('AllKeepable', SuperLocusAll, TranscriptLevelNice, FinalFeatures)
