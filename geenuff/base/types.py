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


########
# GFF
########

# SuperLocus
CODING_GENE = 'coding_gene'
NON_CODING_GENE = 'non_coding_gene'
PSEUDOGENE = 'pseudogene'
OPERON = 'operon'
NCRNA_GENE = 'ncRNA_gene'
BIDIRECTIONAL_PROMOTER_LNCRNA = 'bidirectional_promoter_lncRNA'
TRANSPOSABLE_ELEMENT_GENE = 'transposable_element_gene'
GENE = 'gene'

SuperLocus = make_enum('SuperLocus', CODING_GENE, NON_CODING_GENE, PSEUDOGENE, OPERON, NCRNA_GENE,
                       BIDIRECTIONAL_PROMOTER_LNCRNA, TRANSPOSABLE_ELEMENT_GENE)
SuperLocusHistorical = make_enum('SuperLocusHistorical', GENE)
SuperLocusAll = join_to_enum('SuperLocusAll', SuperLocus, SuperLocusHistorical)

# Transcript
LNC_RNA = 'lnc_RNA'
MRNA = 'mRNA'
MIRNA = 'miRNA'
NCRNA = 'ncRNA'
PIRNA = 'piRNA'
PRE_MIRNA = 'pre_miRNA'
RNASE_MRP_RNA = 'RNase_MRP_RNA'
RRNA = 'rRNA'
SCRNA = 'scRNA'
SNORNA = 'snoRNA'
SNRNA = 'snRNA'
SRP_RNA = 'SRP_RNA'
TRNA = 'tRNA'
Y_RNA = 'Y_RNA'
TRANSCRIPT = 'transcript'
PRIMARY_TRANSCRIPT = 'primary_transcript'
PSEUDOGENIC_TRANSCRIPT = 'pseudogenic_transcript'
UNCONFIRMED_TRANSCRIPT = 'unconfirmed_transcript'
VAULTRNA_PRIMARY_TRANSCRIPT = 'vaultRNA_primary_transcript'
THREE_PRIME_OVERLAPPING_NCRNA = 'three_prime_overlapping_ncrna'
C_GENE_SEGMENT = 'C_gene_segment'
V_GENE_SEGMENT = 'V_gene_segment'
D_GENE_SEGMENT = 'D_gene_segment'
J_GENE_SEGMENT = 'J_gene_segment'
GENE_SEGMENT = 'gene_segment'
TRANSPOSABLE_ELEMENT = 'transposable_element'

TranscriptLevel = make_enum('TranscriptLevel', LNC_RNA, MRNA, MIRNA, NCRNA, PIRNA, PRE_MIRNA,
                            RNASE_MRP_RNA, RRNA, SCRNA, SNORNA, SNRNA, SRP_RNA, TRNA, Y_RNA,
                            TRANSCRIPT, PRIMARY_TRANSCRIPT, PSEUDOGENIC_TRANSCRIPT,
                            UNCONFIRMED_TRANSCRIPT, VAULTRNA_PRIMARY_TRANSCRIPT,
                            THREE_PRIME_OVERLAPPING_NCRNA, C_GENE_SEGMENT, V_GENE_SEGMENT,
                            D_GENE_SEGMENT, J_GENE_SEGMENT, GENE_SEGMENT, TRANSPOSABLE_ELEMENT)

# ignorable because of reasons
SUPERCONTIG = 'supercontig'
CONTIG = 'contig'
CHROMOSOME = 'chromosome'
BIOLOGICAL_REGION = 'biological_region'
REGION = 'region'
SCAFFOLD = 'scaffold'
FIVE_PRIME_UTR = 'five_prime_UTR'
THREE_PRIME_UTR = 'three_prime_UTR'
MATCH = 'match'
CDNA_MATCH = 'cDNA_match'

# ignorable from RefSeq genomes
ANTISENSE_RNA = 'antisense_RNA'
CAAT_SIGNAL = 'CAAT_signal'
CAGE_CLUSTER = 'CAGE_cluster'
CENTROMERE = 'centromere'
CHROMOSOME_BREAKPOINT = 'chromosome_breakpoint'
CONSERVED_REGION = 'conserved_region'
DIRECT_REPEAT = 'direct_repeat'
DISPERSED_REPEAT = 'dispersed_repeat'
D_LOOP = 'D_loop'
DNASEI_HYPERSENSITIVE_SITE = 'DNaseI_hypersensitive_site'
DNASEI_HYPERSENSITIVE_SITE2 = 'DNAseI_hypersensitive_site'
ENHANCER = 'enhancer'
ENHANCER_BLOCKING_ELEMENT = 'enhancer_blocking_element'
EPIGENETICALLY_MODIFIED_REGION = 'epigenetically_modified_region'
GC_RICH_PROMOTER_REGION = 'GC_rich_promoter_region'
GUIDE_RNA = 'guide_RNA'
IMPRINTING_CONTROL_REGION = 'imprinting_control_region'
INSULATOR = 'insulator'
LOCUS_CONTROL_REGION = 'locus_control_region'
MATRIX_ATTACHMENT_SITE = 'matrix_attachment_site'
MEIOTIC_RECOMBINATION_REGION = 'meiotic_recombination_region'
MICROSATELLITE = 'microsatellite'
MINISATELLITE = 'minisatellite'
MITOTIC_RECOMBINATION_REGION = 'mitotic_recombination_region'
MOBILE_GENETIC_ELEMENT = 'mobile_genetic_element'
NON_ALLELIC_HOMOLOGOUS_RECOMBINATION_REGION = 'non_allelic_homologous_recombination_region'
NUCLEOTIDE_CLEAVAGE_SITE = 'nucleotide_cleavage_site'
NUCLEOTIDE_MOTIF = 'nucleotide_motif'
ORIGIN_OF_REPLICATION = 'origin_of_replication'
PROMOTER = 'promoter'
PROTEIN_BINDING_SITE = 'protein_binding_site'
RECOMBINATION_FEATURE = 'recombination_feature'
REGULATORY_REGION = 'regulatory_region'
REPEAT_INSTABILITY_REGION = 'repeat_instability_region'
REPEAT_REGION = 'repeat_region'
REPLICATION_REGULATORY_REGION = 'replication_regulatory_region'
REPLICATION_START_SITE = 'replication_start_site'
RESPONSE_ELEMENT = 'response_element'
RNASE_P_RNA = 'RNase_P_RNA'
SEQUENCE_ALTERATION = 'sequence_alteration'
SEQUENCE_COMPARISON = 'sequence_comparison'
SEQUENCE_FEATURE = 'sequence_feature'
SEQUENCE_SECONDARY_STRUCTURE = 'sequence_secondary_structure'
SILENCER = 'silencer'
TANDEM_REPEAT = 'tandem_repeat'
TATA_box = 'TATA_box'
TELOMERASE_RNA = 'telomerase_RNA'
TRANSCRIPTIONAL_CIS_REGULATORY_REGION = 'transcriptional_cis_regulatory_region'
TSS = 'TSS'
VAULT_RNA = 'vault_RNA'

# ignorable from Augustus
START_CODON = 'start_codon'
STOP_CODON = 'stop_codon'
INTRON = 'intron'
TRANSCRIPTION_START_SITE = 'transcription_start_site'  # transcription_start_site
TRANSCRIPTION_TERMINATION_SITE = 'transcription_end_site'  # transcription_termination_site
TRANSCRIPTION_START_SITE2 = 'tss'  # transcription_start_site (older Augustus version)
TRANSCRIPTION_TERMINATION_SITE2 = 'tts'  # transcription_termination_site
FIVE_PRIME_UTR_LOWER = 'five_prime_utr'
THREE_PRIME_UTR_LOWER = 'three_prime_utr'

IgnorableGFFFeatures = make_enum('IgnorableGFFFeatures', SUPERCONTIG, CONTIG, CHROMOSOME,
                                 BIOLOGICAL_REGION, REGION, SCAFFOLD, FIVE_PRIME_UTR, THREE_PRIME_UTR,
                                 MATCH, CDNA_MATCH, ANTISENSE_RNA, CAAT_SIGNAL, CAGE_CLUSTER, CENTROMERE,
                                 CHROMOSOME_BREAKPOINT, CONSERVED_REGION, DIRECT_REPEAT, DISPERSED_REPEAT,
                                 D_LOOP, DNASEI_HYPERSENSITIVE_SITE, DNASEI_HYPERSENSITIVE_SITE2, ENHANCER,
                                 ENHANCER_BLOCKING_ELEMENT, EPIGENETICALLY_MODIFIED_REGION, GC_RICH_PROMOTER_REGION,
                                 GUIDE_RNA, IMPRINTING_CONTROL_REGION, INSULATOR, LOCUS_CONTROL_REGION,
                                 MATRIX_ATTACHMENT_SITE, MEIOTIC_RECOMBINATION_REGION, MICROSATELLITE, MINISATELLITE,
                                 MITOTIC_RECOMBINATION_REGION, MOBILE_GENETIC_ELEMENT,
                                 NON_ALLELIC_HOMOLOGOUS_RECOMBINATION_REGION, NUCLEOTIDE_CLEAVAGE_SITE, NUCLEOTIDE_MOTIF,
                                 ORIGIN_OF_REPLICATION, PROMOTER, PROTEIN_BINDING_SITE, RECOMBINATION_FEATURE,
                                 REGULATORY_REGION, REPEAT_INSTABILITY_REGION, REPEAT_REGION, REPLICATION_REGULATORY_REGION,
                                 REPLICATION_START_SITE, RESPONSE_ELEMENT, RNASE_P_RNA, SEQUENCE_ALTERATION,
                                 SEQUENCE_COMPARISON, SEQUENCE_FEATURE, SEQUENCE_SECONDARY_STRUCTURE, SILENCER,
                                 TANDEM_REPEAT, TATA_box, TELOMERASE_RNA, TRANSCRIPTIONAL_CIS_REGULATORY_REGION,
                                 TSS, VAULT_RNA, START_CODON, STOP_CODON, INTRON, TRANSCRIPTION_START_SITE,
                                 TRANSCRIPTION_TERMINATION_SITE, TRANSCRIPTION_START_SITE2,
                                 TRANSCRIPTION_TERMINATION_SITE2, FIVE_PRIME_UTR_LOWER,
                                 THREE_PRIME_UTR_LOWER)

# other useful features
EXON = 'exon'
CDS = 'CDS'

UsefulGFFSequenceFeatures = make_enum('UsefulGFFSequenceFeatures', EXON, CDS)
UsefulGFFFeatures = join_to_enum('UsefulGFFFeatures', SuperLocusAll, TranscriptLevel,
                                 UsefulGFFSequenceFeatures)
AllKnownGFFFeatures = join_to_enum('AllKnownGFFFeatures', IgnorableGFFFeatures, UsefulGFFFeatures)

########
# Geenuff
########

GEENUFF_TRANSCRIPT = 'geenuff_transcript'
GEENUFF_CDS= 'geenuff_cds'
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
