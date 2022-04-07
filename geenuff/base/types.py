from .helpers import (join_to_enum,
                      make_enum,
                      join_to_enum_strip_redundancy,
                      subtract_enum)


from .so import (SOSequenceFeatures,
                 SOSuperLocusFeatures,
                 SOTranscriptFeatures,
                 SOCDSFeatures,
                 SOExonFeatures)
########
# GFF
########

# SuperLocus
SuperLocusAll = join_to_enum('SuperLocusAll', SOSuperLocusFeatures)

# Transcript
TranscriptLevel = join_to_enum('TranscriptLevel', SOTranscriptFeatures)

# ignorable (and non-SO sequence features) from Augustus
TRANSCRIPTION_START_SITE = 'transcription_start_site'  # transcription_start_site
TRANSCRIPTION_START_SITE2 = 'tss'  # transcription_start_site (older Augustus version)
TRANSCRIPTION_TERMINATION_SITE2 = 'tts'  # transcription_termination_site
FIVE_PRIME_UTR_LOWER = 'five_prime_utr'
THREE_PRIME_UTR_LOWER = 'three_prime_utr'

IgnorableGFFFeatures = make_enum('IgnorableGFFFeatures', TRANSCRIPTION_START_SITE,
                                 TRANSCRIPTION_START_SITE2,
                                 TRANSCRIPTION_TERMINATION_SITE2, FIVE_PRIME_UTR_LOWER,
                                 THREE_PRIME_UTR_LOWER)

# other useful features
ExonLevel = join_to_enum('ExonLevel', SOExonFeatures)
CDSLevel = join_to_enum('CDSLevel', SOCDSFeatures)

UsefulGFFFeatures = join_to_enum('UsefulGFFFeatures', SuperLocusAll, TranscriptLevel,
                                 ExonLevel, CDSLevel)


# we will need a list that can be ignored, so w/o things like exon and CDS, but that we can nevertheless
# use to check for typos/blatantly invalid gff3
OnlySoFeatures = subtract_enum('OnlySoFeatures', subtract_from=SOSequenceFeatures, subtract_this=UsefulGFFFeatures)
IgnorableGFFFeatures = join_to_enum('IgnorableGffFeatures', IgnorableGFFFeatures, OnlySoFeatures)

AllKnownGFFFeatures = join_to_enum('AllKnownGFFFeatures', IgnorableGFFFeatures, UsefulGFFFeatures)

###########
# GeenuFF #
###########

GEENUFF_TRANSCRIPT = 'geenuff_transcript'
GEENUFF_CDS = 'geenuff_cds'
GEENUFF_INTRON = 'geenuff_intron'
GeenuffSequenceFeature = make_enum('GeenuffSequenceFeature', GEENUFF_TRANSCRIPT, GEENUFF_CDS,
                                   GEENUFF_INTRON)

# Geenuff error types
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
MISMATCHING_STRANDS = 'missmatching_strands'
TRUNCATED_INTRON = 'truncated_intron'
Errors = make_enum('Errors', MISSING_UTR_5P, MISSING_UTR_3P, EMPTY_SUPER_LOCUS, MISSING_START_CODON,
                   MISSING_STOP_CODON, WRONG_PHASE_5P, MISMATCHED_PHASE_3P, OVERLAPPING_EXONS,
                   TOO_SHORT_INTRON, SL_OVERLAP_ERROR, MISMATCHING_STRANDS, TRUNCATED_INTRON)

# frequently used for lookups
geenuff_error_type_values = [t.value for t in Errors]

GeenuffFeature = join_to_enum('GeenuffFeature', GeenuffSequenceFeature, Errors)
