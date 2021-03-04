import os
import pytest
from ..applications.importer import ImportController
from ..applications.exporters.sequence import FastaExportController
from ..applications.exporters.lengths import LengthExportController
from ..applications.exporters.json import JsonExportController, FeatureJsonable, TranscriptJsonable, SuperLocusJsonable
from ..applications.exporter import MODES
from geenuff.base import orm, types
import json
from geenuff.applications.exporter import GeenuffExportController


EXPORTING_PFX = 'testdata/exporting.'
EXONEXONCDS_PFX = 'testdata/exonexonCDS.'

EXPORTING_DB = EXPORTING_PFX + 'sqlite3'
EXONEXONCDS_DB = EXONEXONCDS_PFX + 'sqlite3'
SECOND_EXP_DB = 'testdata/tmp_ex2.sqlite3'


@pytest.fixture(scope="module", autouse=True)
def prepare_and_cleanup():
    if not os.getcwd().endswith('GeenuFF/geenuff'):
        pytest.exit('Tests need to be run from GeenuFF/geenuff directory')

    for pfx in [EXPORTING_PFX, EXONEXONCDS_PFX]:
        if not os.path.exists(pfx + 'sqlite3'):
            controller = ImportController(database_path='sqlite:///' + pfx + 'sqlite3')
            controller.add_genome(pfx + 'fa', pfx + 'gff3', clean_gff=True,
                                  genome_args={'species': 'dummy'})

    # add three, so we have some to include / exclude for 'test_exporter_genomes_or_exclude'
    if not os.path.exists(SECOND_EXP_DB):
        controller = ImportController(database_path='sqlite:///' + SECOND_EXP_DB)
        controller.add_genome('testdata/exporting.fa', 'testdata/exporting.gff3', clean_gff=True,
                              genome_args={'species': 'exporting'})
        controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff',
                              genome_args={'species': 'dummyloci'})
        controller.add_genome('testdata/exonexonCDS.fa', 'testdata/exonexonCDS.gff3',
                              genome_args={'species': 'exonexonCDS'})

    yield
    for db in [EXPORTING_DB, EXONEXONCDS_DB]:
        os.remove(db)

    os.remove(SECOND_EXP_DB)


def seq_len_controllers(mode, longest=False, db=EXPORTING_DB):
    econtroller = FastaExportController(db_path_in='sqlite:///' + db, longest=longest)
    econtroller.prep_ranges(range_function=MODES[mode])

    lcontroller = LengthExportController(db_path_in='sqlite:///' + db, longest=longest)
    lcontroller.prep_ranges(range_function=MODES[mode])
    return econtroller, lcontroller


def compare2controllers(expect, econtroller, lcontroller):
    iexpect = iter(expect)
    ilengths = iter(lcontroller.export_ranges)
    for i, grp in enumerate(econtroller.export_ranges):
        exp = next(iexpect)
        lgrp = next(ilengths)
        seq = econtroller.get_seq(grp)
        length = lcontroller.get_length(lgrp)
        try:
            assert ''.join(seq) == exp[2]
            assert length == exp[1]
        except AssertionError as e:
            print(exp[0])
            print(i)
            raise e


def test_get_intron_seqs():
    """checks that expected intron sequences are produced"""
    expect = [("Chr1:195000-199000:1543-2114", 572,
               "GTACTTCCAAATCTTCAATTTTGATTCTAAAGATTGGTCCTTTTACTCTGTTTCTCAATT"
               "TGAGTTTTAGGTATTCTTTGATTTTGTATTGGTTTCATTCTAAATATTCATCCTTTACTC"
               "AACTTCTAGATAAGGGATTTAGGTATTCTCAAATTTCCGATTTGATTCCTTTACTCGTTT"
               "CTAGATTGGGGTTTTAGGAATTACCAGTTGGGGGTTTTGCAATTTGCGTAATCAAAGAAT"
               "TTTATTTGTTGTATTGCTTGGTATTGAAGTTTGTCTCTGTTTCTCTACCTCGTCATGTAA"
               "TGTGCTTAGATCCATTAAGTAAATGCTTGTGGATATTTATGTAGATGGTTAAGAGTGATC"
               "GTGATCAGAGTCCTTCTCTTATTTAACTGCATTGCCTGTGAGTTGTGGTCCTGAAGGTTG"
               "TTGTTATTATTGAATTCTATGTATGTATAGATTATGTCATTGGTCTCATGTGGTTTTTAT"
               "GGGTAACGTCTTTACTAATAATAGCACTATGCTTCTGGATTTTGATCTATGTGATCTGTA"
               "ACATTTCTAGTTGGTGTGTCTTTGATTGCCAG"),
              ("Chr1:195000-199000:2220-2299", 80,
               "GTATATATACCGCTGCTCGTATCTCTTTTCCGGTGTTACAAAAGCGATGTCGTGACCTAA"
               "TGCTGGGTTCGTTACTATAG"),
              ("Chr1:195000-199000:2423-2500", 78,
               "GTAAGTCTGGAATAGCTTTTGAGTTGTCCTCTATGTTTATAAGCTATTGTTGTGTGTAAA"
               "CCTTTGTTATATCTGTAG"),
              ("Chr1:195000-199000:2672-2774", 103,
               "GTAAACTATTAAACTCATTAACTCTCTCCTGCAATCTGCAAGGCAGTCTTTAGGAATGTG"
               "AATATTAGGAAATAACTTTTACTTTGTGGGTTGATTTGTTTAG"),
              ("Chr1:195000-199000:2905-2973", 69,
               "GTATTTGATAAGCATATTCAATCCTCTCTATGTTTTTGTAAATGGTCTTACTATGTGTCG"
               "TGGCTGTAG"),
              ("Chr1:195000-199000:3184-3266", 83,
               "GTATACAAAACAATTTGCCTTTACGTTTTTACATTTCTTTAAGAGTTTGAAACATGTCTA"
               "AAGCTGGGATAATATTTTTGCAG")]
    expect += expect[:4] + expect[5:6]

    # test data has been simplified from an augustus run that previously resulted in erroneous masks
    # and from a partial gene model in the Rcommumnis genome
    econtroller, lcontroller = seq_len_controllers('introns')
    # expect [(fa_id (samtools faidx), length, sequence)]
    # and the fa_id is just for knowing where it went wrong if it does
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 11
    compare2controllers(expect, econtroller, lcontroller)
    # and now just the longest (1st & 3rd)
    econtroller, lcontroller = seq_len_controllers('introns', longest=True)
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 6
    compare2controllers(expect[:6], econtroller, lcontroller)


def test_get_exons():
    """checks that individual expected exon sequences are produced"""
    expect = [("Chr1:195000-199000:780-1542", 763,
               "AAGCCTTTCTCTTTAAATTCGTTATCGTTTTTTTTATTTTATCAATTTAATCTTTTTATT"
               "TGTTTTGTTCTTCCTCGATTCAACACTCGATGCTGTGACAAAGCCCAGATTTTGGCCGGA"
               "AAGTTTTGTGTGTTTCCGGCGAAGAATGCGAAGAAACCGGAAGAATCTGTAACGGAATCT"
               "AAGCTAAAAGTTAAAGTACGATGTCGGTGTACGACGCTGCTTTCCTTAATACAGAGCTTT"
               "CGAAACCGACATCGATTTTTGGTCTCCGGCTATGGGTCGTGATCGGAATCTTACTTGGAT"
               "CTCTAATTGTCATCGCACTCTTTCTTCTCTCCCTCTGCTTAACTTCTCGCCGGAAAAATC"
               "GAAAGCCGAGAGCCGATTTCGCCTCCGCCGCCATCGCTACACCGCCGATTTCAAAGGAGA"
               "TTAAAGAGATCGTTCCGGCGCAGAATCAGTCTGTTCCGGCGGAGATCCAGGTCGATATCG"
               "GGAAGATCGAGCATCGAGTGGTGTTTTCAGATCGAGTGTCGAGTGGTGAGAGTAGAGGAA"
               "CAGCGAGTGCAAGTGAAACGGCGTCGTATTCCGGTAGCGGGAATTGTGGGCCGGAGGTGT"
               "CGCATCTTGGATGGGGCCGATGGTATACTCTGAGAGAGCTTGAAGCGGCCACGAATGGGC"
               "TTTGTGAAGAGAATGTAATCGGAGAAGGTGGTTACGGGATTGTGTATCGTGGCATTTTAA"
               "CTGATGGAACCAAAGTCGCCGTCAAGAACTTGCTTAACAATAG"),
              ("Chr1:195000-199000:2115-2219", 105,
               "GGGTCAAGCAGAGAAGGAATTCAAAGTAGAAGTGGAAGTCATTGGGCGAGTACGACACAA"
               "GAATCTTGTTAGGCTTTTAGGGTATTGCGTGGAAGGTGCATACAG"),
              ("Chr1:195000-199000:2300-2422", 123,
               "GATGCTCGTGTATGACTTTGTCGACAATGGTAATTTGGAGCAATGGATTCACGGTGATGT"
               "TGGCGATGTCAGCCCGCTAACTTGGGATATACGTATGAATATTATACTGGGGATGGCCAA"
               "AGG"),
              ("Chr1:195000-199000:2501-2671", 171,
               "ATTGGCGTATCTACACGAGGGTCTTGAGCCAAAAGTTGTTCATCGGGATATTAAATCAAG"
               "CAATATCTTACTTGATCGCCAATGGAATGCTAAGGTTTCGGATTTTGGACTTGCTAAGCT"
               "CTTGGGGTCTGAGAGCAGTTATGTGACTACTCGTGTGATGGGAACTTTCGG"),
              ("Chr1:195000-199000:2775-2904", 130,
               "TTATGTAGCACCAGAATACGCTTGCACCGGAATGTTAAACGAGAAGAGTGATATCTATAG"
               "CTTCGGAATACTAATCATGGAGATAATCACTGGAAGAAACCCGGTTGATTATAGTCGGCC"
               "TCAAGGAGAG"),
              ("Chr1:195000-199000:2974-3183", 210,
               "ACAAATCTAGTGGATTGGCTTAAATCAATGGTGGGAAACCGAAGATCAGAAGAAGTTGTT"
               "GATCCGAAAATACCCGAACCACCATCCTCAAAGGCTCTTAAACGAGTGTTGCTTGTAGCT"
               "TTGCGTTGTGTGGATCCTGATGCGAACAAGAGACCTAAAATGGGTCATATCATACATATG"
               "CTTGAAGCCGAAGATCTACTCTATCGCGAT"),
              ("Chr1:195000-199000:3267-3684", 418,
               "GAACGCCGAACAACAAGGGACCATGGAAGCCGCGAGAGACAAGAGACAGCTGTGGTTGCT"
               "GCCGGTAGTGAAAGTGGTGAGAGCGGTTCACGGCATCATCAGCAAAAGCAAAGATGAAAA"
               "AAGAGAGTCACTTGGGTTAAGTGATTTCCACACGACCATTATATTTCATTTTTTTTTCCG"
               "GTATATTATTGTCTCACTTACATAATATTTTCATTCTTTGGCTCTTCGAGTCTTGGCCTT"
               "ACGAAATTAACATTCTTCTTCTTTGTAACCTTTTGTCAACGATCACATTCACATGTCACT"
               "TATTAGATATTGTAATATGTAATGTTTGGACCGACGTCGAAGCATAAGCAGGTGATTGGT"
               "CGATGGATCAAACTTATCGCTTTGGACGAAATGGACCACAATGATTCTTTTTTAGCTC"),
              ("Chr1:195000-199000:812-1542", 731,
               "TTTATTTTATCAATTTAATCTTTTTATTTGTTTTGTTCTTCCTCGATTCAACACTCGATG"
               "CTGTGACAAAGCCCAGATTTTGGCCGGAAAGTTTTGTGTGTTTCCGGCGAAGAATGCGAA"
               "GAAACCGGAAGAATCTGTAACGGAATCTAAGCTAAAAGTTAAAGTACGATGTCGGTGTAC"
               "GACGCTGCTTTCCTTAATACAGAGCTTTCGAAACCGACATCGATTTTTGGTCTCCGGCTA"
               "TGGGTCGTGATCGGAATCTTACTTGGATCTCTAATTGTCATCGCACTCTTTCTTCTCTCC"
               "CTCTGCTTAACTTCTCGCCGGAAAAATCGAAAGCCGAGAGCCGATTTCGCCTCCGCCGCC"
               "ATCGCTACACCGCCGATTTCAAAGGAGATTAAAGAGATCGTTCCGGCGCAGAATCAGTCT"
               "GTTCCGGCGGAGATCCAGGTCGATATCGGGAAGATCGAGCATCGAGTGGTGTTTTCAGAT"
               "CGAGTGTCGAGTGGTGAGAGTAGAGGAACAGCGAGTGCAAGTGAAACGGCGTCGTATTCC"
               "GGTAGCGGGAATTGTGGGCCGGAGGTGTCGCATCTTGGATGGGGCCGATGGTATACTCTG"
               "AGAGAGCTTGAAGCGGCCACGAATGGGCTTTGTGAAGAGAATGTAATCGGAGAAGGTGGT"
               "TACGGGATTGTGTATCGTGGCATTTTAACTGATGGAACCAAAGTCGCCGTCAAGAACTTG"
               "CTTAACAATAG")]
    expect += expect[1:4]
    expect += [("Chr1:195000-199000:2775-3183", 409,
                "TTATGTAGCACCAGAATACGCTTGCACCGGAATGTTAAACGAGAAGAGTGATATCTATAG"
                "CTTCGGAATACTAATCATGGAGATAATCACTGGAAGAAACCCGGTTGATTATAGTCGGCC"
                "TCAAGGAGAGGTATTTGATAAGCATATTCAATCCTCTCTATGTTTTTGTAAATGGTCTTA"
                "CTATGTGTCGTGGCTGTAGACAAATCTAGTGGATTGGCTTAAATCAATGGTGGGAAACCG"
                "AAGATCAGAAGAAGTTGTTGATCCGAAAATACCCGAACCACCATCCTCAAAGGCTCTTAA"
                "ACGAGTGTTGCTTGTAGCTTTGCGTTGTGTGGATCCTGATGCGAACAAGAGACCTAAAAT"
                "GGGTCATATCATACATATGCTTGAAGCCGAAGATCTACTCTATCGCGAT"),
               ("Chr1:195000-199000:3267-3635", 369,
                "GAACGCCGAACAACAAGGGACCATGGAAGCCGCGAGAGACAAGAGACAGCTGTGGTTGCT"
                "GCCGGTAGTGAAAGTGGTGAGAGCGGTTCACGGCATCATCAGCAAAAGCAAAGATGAAAA"
                "AAGAGAGTCACTTGGGTTAAGTGATTTCCACACGACCATTATATTTCATTTTTTTTTCCG"
                "GTATATTATTGTCTCACTTACATAATATTTTCATTCTTTGGCTCTTCGAGTCTTGGCCTT"
                "ACGAAATTAACATTCTTCTTCTTTGTAACCTTTTGTCAACGATCACATTCACATGTCACT"
                "TATTAGATATTGTAATATGTAATGTTTGGACCGACGTCGAAGCATAAGCAGGTGATTGGT"
                "CGATGGATC"),
               ("27488:1-3000:2777-2962", 186,
                "TGGTTATGAAAATAATAAGAAATTTAAAGTTAAATGTCTACATGAACATGTTGAAAATCA"
                "CCTTAAGACTGTAAAAAATGAACGAGCAACTATTTCAACACTTCAAGCTAAAAATGGATT"
                "TAGATGGGATGATACTGCAACAAAATATGTTTATGATGAACAAATGATGGTAAGCATTCT"
                "TAATAA")]

    econtroller, lcontroller = seq_len_controllers('exons')
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 14
    compare2controllers(expect, econtroller, lcontroller)
    # and now just the longest (1st & 3rd)
    econtroller, lcontroller = seq_len_controllers('exons', longest=True)
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 8
    expect = expect[:7] + expect[-1:]
    compare2controllers(expect, econtroller, lcontroller)


def test_get_mRNA():
    """tests that exons are propperly spliced into mRNA sequence"""
    expect = [("AT1G01540.2.TAIR10", 1920,
               "AAGCCTTTCTCTTTAAATTCGTTATCGTTTTTTTTATTTTATCAATTTAATCTTTTTATTTGTTTTGTTC"
               "TTCCTCGATTCAACACTCGATGCTGTGACAAAGCCCAGATTTTGGCCGGAAAGTTTTGTGTGTTTCCGGC"
               "GAAGAATGCGAAGAAACCGGAAGAATCTGTAACGGAATCTAAGCTAAAAGTTAAAGTACGATGTCGGTGT"
               "ACGACGCTGCTTTCCTTAATACAGAGCTTTCGAAACCGACATCGATTTTTGGTCTCCGGCTATGGGTCGT"
               "GATCGGAATCTTACTTGGATCTCTAATTGTCATCGCACTCTTTCTTCTCTCCCTCTGCTTAACTTCTCGC"
               "CGGAAAAATCGAAAGCCGAGAGCCGATTTCGCCTCCGCCGCCATCGCTACACCGCCGATTTCAAAGGAGA"
               "TTAAAGAGATCGTTCCGGCGCAGAATCAGTCTGTTCCGGCGGAGATCCAGGTCGATATCGGGAAGATCGA"
               "GCATCGAGTGGTGTTTTCAGATCGAGTGTCGAGTGGTGAGAGTAGAGGAACAGCGAGTGCAAGTGAAACG"
               "GCGTCGTATTCCGGTAGCGGGAATTGTGGGCCGGAGGTGTCGCATCTTGGATGGGGCCGATGGTATACTC"
               "TGAGAGAGCTTGAAGCGGCCACGAATGGGCTTTGTGAAGAGAATGTAATCGGAGAAGGTGGTTACGGGAT"
               "TGTGTATCGTGGCATTTTAACTGATGGAACCAAAGTCGCCGTCAAGAACTTGCTTAACAATAGGGGTCAA"
               "GCAGAGAAGGAATTCAAAGTAGAAGTGGAAGTCATTGGGCGAGTACGACACAAGAATCTTGTTAGGCTTT"
               "TAGGGTATTGCGTGGAAGGTGCATACAGGATGCTCGTGTATGACTTTGTCGACAATGGTAATTTGGAGCA"
               "ATGGATTCACGGTGATGTTGGCGATGTCAGCCCGCTAACTTGGGATATACGTATGAATATTATACTGGGG"
               "ATGGCCAAAGGATTGGCGTATCTACACGAGGGTCTTGAGCCAAAAGTTGTTCATCGGGATATTAAATCAA"
               "GCAATATCTTACTTGATCGCCAATGGAATGCTAAGGTTTCGGATTTTGGACTTGCTAAGCTCTTGGGGTC"
               "TGAGAGCAGTTATGTGACTACTCGTGTGATGGGAACTTTCGGTTATGTAGCACCAGAATACGCTTGCACC"
               "GGAATGTTAAACGAGAAGAGTGATATCTATAGCTTCGGAATACTAATCATGGAGATAATCACTGGAAGAA"
               "ACCCGGTTGATTATAGTCGGCCTCAAGGAGAGACAAATCTAGTGGATTGGCTTAAATCAATGGTGGGAAA"
               "CCGAAGATCAGAAGAAGTTGTTGATCCGAAAATACCCGAACCACCATCCTCAAAGGCTCTTAAACGAGTG"
               "TTGCTTGTAGCTTTGCGTTGTGTGGATCCTGATGCGAACAAGAGACCTAAAATGGGTCATATCATACATA"
               "TGCTTGAAGCCGAAGATCTACTCTATCGCGATGAACGCCGAACAACAAGGGACCATGGAAGCCGCGAGAG"
               "ACAAGAGACAGCTGTGGTTGCTGCCGGTAGTGAAAGTGGTGAGAGCGGTTCACGGCATCATCAGCAAAAG"
               "CAAAGATGAAAAAAGAGAGTCACTTGGGTTAAGTGATTTCCACACGACCATTATATTTCATTTTTTTTTC"
               "CGGTATATTATTGTCTCACTTACATAATATTTTCATTCTTTGGCTCTTCGAGTCTTGGCCTTACGAAATT"
               "AACATTCTTCTTCTTTGTAACCTTTTGTCAACGATCACATTCACATGTCACTTATTAGATATTGTAATAT"
               "GTAATGTTTGGACCGACGTCGAAGCATAAGCAGGTGATTGGTCGATGGATCAAACTTATCGCTTTGGACG"
               "AAATGGACCACAATGATTCTTTTTTAGCTC"),
              ("AT1G01540.1.TAIR10", 1908,
               "TTTATTTTATCAATTTAATCTTTTTATTTGTTTTGTTCTTCCTCGATTCAACACTCGATGCTGTGACAAA"
               "GCCCAGATTTTGGCCGGAAAGTTTTGTGTGTTTCCGGCGAAGAATGCGAAGAAACCGGAAGAATCTGTAA"
               "CGGAATCTAAGCTAAAAGTTAAAGTACGATGTCGGTGTACGACGCTGCTTTCCTTAATACAGAGCTTTCG"
               "AAACCGACATCGATTTTTGGTCTCCGGCTATGGGTCGTGATCGGAATCTTACTTGGATCTCTAATTGTCA"
               "TCGCACTCTTTCTTCTCTCCCTCTGCTTAACTTCTCGCCGGAAAAATCGAAAGCCGAGAGCCGATTTCGC"
               "CTCCGCCGCCATCGCTACACCGCCGATTTCAAAGGAGATTAAAGAGATCGTTCCGGCGCAGAATCAGTCT"
               "GTTCCGGCGGAGATCCAGGTCGATATCGGGAAGATCGAGCATCGAGTGGTGTTTTCAGATCGAGTGTCGA"
               "GTGGTGAGAGTAGAGGAACAGCGAGTGCAAGTGAAACGGCGTCGTATTCCGGTAGCGGGAATTGTGGGCC"
               "GGAGGTGTCGCATCTTGGATGGGGCCGATGGTATACTCTGAGAGAGCTTGAAGCGGCCACGAATGGGCTT"
               "TGTGAAGAGAATGTAATCGGAGAAGGTGGTTACGGGATTGTGTATCGTGGCATTTTAACTGATGGAACCA"
               "AAGTCGCCGTCAAGAACTTGCTTAACAATAGGGGTCAAGCAGAGAAGGAATTCAAAGTAGAAGTGGAAGT"
               "CATTGGGCGAGTACGACACAAGAATCTTGTTAGGCTTTTAGGGTATTGCGTGGAAGGTGCATACAGGATG"
               "CTCGTGTATGACTTTGTCGACAATGGTAATTTGGAGCAATGGATTCACGGTGATGTTGGCGATGTCAGCC"
               "CGCTAACTTGGGATATACGTATGAATATTATACTGGGGATGGCCAAAGGATTGGCGTATCTACACGAGGG"
               "TCTTGAGCCAAAAGTTGTTCATCGGGATATTAAATCAAGCAATATCTTACTTGATCGCCAATGGAATGCT"
               "AAGGTTTCGGATTTTGGACTTGCTAAGCTCTTGGGGTCTGAGAGCAGTTATGTGACTACTCGTGTGATGG"
               "GAACTTTCGGTTATGTAGCACCAGAATACGCTTGCACCGGAATGTTAAACGAGAAGAGTGATATCTATAG"
               "CTTCGGAATACTAATCATGGAGATAATCACTGGAAGAAACCCGGTTGATTATAGTCGGCCTCAAGGAGAG"
               "GTATTTGATAAGCATATTCAATCCTCTCTATGTTTTTGTAAATGGTCTTACTATGTGTCGTGGCTGTAGA"
               "CAAATCTAGTGGATTGGCTTAAATCAATGGTGGGAAACCGAAGATCAGAAGAAGTTGTTGATCCGAAAAT"
               "ACCCGAACCACCATCCTCAAAGGCTCTTAAACGAGTGTTGCTTGTAGCTTTGCGTTGTGTGGATCCTGAT"
               "GCGAACAAGAGACCTAAAATGGGTCATATCATACATATGCTTGAAGCCGAAGATCTACTCTATCGCGATG"
               "AACGCCGAACAACAAGGGACCATGGAAGCCGCGAGAGACAAGAGACAGCTGTGGTTGCTGCCGGTAGTGA"
               "AAGTGGTGAGAGCGGTTCACGGCATCATCAGCAAAAGCAAAGATGAAAAAAGAGAGTCACTTGGGTTAAG"
               "TGATTTCCACACGACCATTATATTTCATTTTTTTTTCCGGTATATTATTGTCTCACTTACATAATATTTT"
               "CATTCTTTGGCTCTTCGAGTCTTGGCCTTACGAAATTAACATTCTTCTTCTTTGTAACCTTTTGTCAACG"
               "ATCACATTCACATGTCACTTATTAGATATTGTAATATGTAATGTTTGGACCGACGTCGAAGCATAAGCAG"
               "GTGATTGGTCGATGGATC"),
              ("27488.m000034.TIGRR0.1", 186,
               "TGGTTATGAAAATAATAAGAAATTTAAAGTTAAATGTCTACATGAACATGTTGAAAATCACCTTAAGACT"
               "GTAAAAAATGAACGAGCAACTATTTCAACACTTCAAGCTAAAAATGGATTTAGATGGGATGATACTGCAA"
               "CAAAATATGTTTATGATGAACAAATGATGGTAAGCATTCTTAATAA")]

    econtroller, lcontroller = seq_len_controllers('mRNA')
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 3
    compare2controllers(expect, econtroller, lcontroller)
    # and now just the longest (1st & 3rd)
    econtroller, lcontroller = seq_len_controllers('mRNA', longest=True)
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 2
    expect = expect[:1] + expect[-1:]
    compare2controllers(expect, econtroller, lcontroller)


def test_get_pre_mRNA():
    """tests returning sequence of raw transcript (transcription start-end)"""
    expect = [("Chr1:195000-199000:780-3684", 2905,
               "AAGCCTTTCTCTTTAAATTCGTTATCGTTTTTTTTATTTTATCAATTTAATCTTTTTATT"
               "TGTTTTGTTCTTCCTCGATTCAACACTCGATGCTGTGACAAAGCCCAGATTTTGGCCGGA"
               "AAGTTTTGTGTGTTTCCGGCGAAGAATGCGAAGAAACCGGAAGAATCTGTAACGGAATCT"
               "AAGCTAAAAGTTAAAGTACGATGTCGGTGTACGACGCTGCTTTCCTTAATACAGAGCTTT"
               "CGAAACCGACATCGATTTTTGGTCTCCGGCTATGGGTCGTGATCGGAATCTTACTTGGAT"
               "CTCTAATTGTCATCGCACTCTTTCTTCTCTCCCTCTGCTTAACTTCTCGCCGGAAAAATC"
               "GAAAGCCGAGAGCCGATTTCGCCTCCGCCGCCATCGCTACACCGCCGATTTCAAAGGAGA"
               "TTAAAGAGATCGTTCCGGCGCAGAATCAGTCTGTTCCGGCGGAGATCCAGGTCGATATCG"
               "GGAAGATCGAGCATCGAGTGGTGTTTTCAGATCGAGTGTCGAGTGGTGAGAGTAGAGGAA"
               "CAGCGAGTGCAAGTGAAACGGCGTCGTATTCCGGTAGCGGGAATTGTGGGCCGGAGGTGT"
               "CGCATCTTGGATGGGGCCGATGGTATACTCTGAGAGAGCTTGAAGCGGCCACGAATGGGC"
               "TTTGTGAAGAGAATGTAATCGGAGAAGGTGGTTACGGGATTGTGTATCGTGGCATTTTAA"
               "CTGATGGAACCAAAGTCGCCGTCAAGAACTTGCTTAACAATAGGTACTTCCAAATCTTCA"
               "ATTTTGATTCTAAAGATTGGTCCTTTTACTCTGTTTCTCAATTTGAGTTTTAGGTATTCT"
               "TTGATTTTGTATTGGTTTCATTCTAAATATTCATCCTTTACTCAACTTCTAGATAAGGGA"
               "TTTAGGTATTCTCAAATTTCCGATTTGATTCCTTTACTCGTTTCTAGATTGGGGTTTTAG"
               "GAATTACCAGTTGGGGGTTTTGCAATTTGCGTAATCAAAGAATTTTATTTGTTGTATTGC"
               "TTGGTATTGAAGTTTGTCTCTGTTTCTCTACCTCGTCATGTAATGTGCTTAGATCCATTA"
               "AGTAAATGCTTGTGGATATTTATGTAGATGGTTAAGAGTGATCGTGATCAGAGTCCTTCT"
               "CTTATTTAACTGCATTGCCTGTGAGTTGTGGTCCTGAAGGTTGTTGTTATTATTGAATTC"
               "TATGTATGTATAGATTATGTCATTGGTCTCATGTGGTTTTTATGGGTAACGTCTTTACTA"
               "ATAATAGCACTATGCTTCTGGATTTTGATCTATGTGATCTGTAACATTTCTAGTTGGTGT"
               "GTCTTTGATTGCCAGGGGTCAAGCAGAGAAGGAATTCAAAGTAGAAGTGGAAGTCATTGG"
               "GCGAGTACGACACAAGAATCTTGTTAGGCTTTTAGGGTATTGCGTGGAAGGTGCATACAG"
               "GTATATATACCGCTGCTCGTATCTCTTTTCCGGTGTTACAAAAGCGATGTCGTGACCTAA"
               "TGCTGGGTTCGTTACTATAGGATGCTCGTGTATGACTTTGTCGACAATGGTAATTTGGAG"
               "CAATGGATTCACGGTGATGTTGGCGATGTCAGCCCGCTAACTTGGGATATACGTATGAAT"
               "ATTATACTGGGGATGGCCAAAGGGTAAGTCTGGAATAGCTTTTGAGTTGTCCTCTATGTT"
               "TATAAGCTATTGTTGTGTGTAAACCTTTGTTATATCTGTAGATTGGCGTATCTACACGAG"
               "GGTCTTGAGCCAAAAGTTGTTCATCGGGATATTAAATCAAGCAATATCTTACTTGATCGC"
               "CAATGGAATGCTAAGGTTTCGGATTTTGGACTTGCTAAGCTCTTGGGGTCTGAGAGCAGT"
               "TATGTGACTACTCGTGTGATGGGAACTTTCGGGTAAACTATTAAACTCATTAACTCTCTC"
               "CTGCAATCTGCAAGGCAGTCTTTAGGAATGTGAATATTAGGAAATAACTTTTACTTTGTG"
               "GGTTGATTTGTTTAGTTATGTAGCACCAGAATACGCTTGCACCGGAATGTTAAACGAGAA"
               "GAGTGATATCTATAGCTTCGGAATACTAATCATGGAGATAATCACTGGAAGAAACCCGGT"
               "TGATTATAGTCGGCCTCAAGGAGAGGTATTTGATAAGCATATTCAATCCTCTCTATGTTT"
               "TTGTAAATGGTCTTACTATGTGTCGTGGCTGTAGACAAATCTAGTGGATTGGCTTAAATC"
               "AATGGTGGGAAACCGAAGATCAGAAGAAGTTGTTGATCCGAAAATACCCGAACCACCATC"
               "CTCAAAGGCTCTTAAACGAGTGTTGCTTGTAGCTTTGCGTTGTGTGGATCCTGATGCGAA"
               "CAAGAGACCTAAAATGGGTCATATCATACATATGCTTGAAGCCGAAGATCTACTCTATCG"
               "CGATGTATACAAAACAATTTGCCTTTACGTTTTTACATTTCTTTAAGAGTTTGAAACATG"
               "TCTAAAGCTGGGATAATATTTTTGCAGGAACGCCGAACAACAAGGGACCATGGAAGCCGC"
               "GAGAGACAAGAGACAGCTGTGGTTGCTGCCGGTAGTGAAAGTGGTGAGAGCGGTTCACGG"
               "CATCATCAGCAAAAGCAAAGATGAAAAAAGAGAGTCACTTGGGTTAAGTGATTTCCACAC"
               "GACCATTATATTTCATTTTTTTTTCCGGTATATTATTGTCTCACTTACATAATATTTTCA"
               "TTCTTTGGCTCTTCGAGTCTTGGCCTTACGAAATTAACATTCTTCTTCTTTGTAACCTTT"
               "TGTCAACGATCACATTCACATGTCACTTATTAGATATTGTAATATGTAATGTTTGGACCG"
               "ACGTCGAAGCATAAGCAGGTGATTGGTCGATGGATCAAACTTATCGCTTTGGACGAAATG"
               "GACCACAATGATTCTTTTTTAGCTC"),
              ("Chr1:195000-199000:812-3635", 2824,
               "TTTATTTTATCAATTTAATCTTTTTATTTGTTTTGTTCTTCCTCGATTCAACACTCGATG"
               "CTGTGACAAAGCCCAGATTTTGGCCGGAAAGTTTTGTGTGTTTCCGGCGAAGAATGCGAA"
               "GAAACCGGAAGAATCTGTAACGGAATCTAAGCTAAAAGTTAAAGTACGATGTCGGTGTAC"
               "GACGCTGCTTTCCTTAATACAGAGCTTTCGAAACCGACATCGATTTTTGGTCTCCGGCTA"
               "TGGGTCGTGATCGGAATCTTACTTGGATCTCTAATTGTCATCGCACTCTTTCTTCTCTCC"
               "CTCTGCTTAACTTCTCGCCGGAAAAATCGAAAGCCGAGAGCCGATTTCGCCTCCGCCGCC"
               "ATCGCTACACCGCCGATTTCAAAGGAGATTAAAGAGATCGTTCCGGCGCAGAATCAGTCT"
               "GTTCCGGCGGAGATCCAGGTCGATATCGGGAAGATCGAGCATCGAGTGGTGTTTTCAGAT"
               "CGAGTGTCGAGTGGTGAGAGTAGAGGAACAGCGAGTGCAAGTGAAACGGCGTCGTATTCC"
               "GGTAGCGGGAATTGTGGGCCGGAGGTGTCGCATCTTGGATGGGGCCGATGGTATACTCTG"
               "AGAGAGCTTGAAGCGGCCACGAATGGGCTTTGTGAAGAGAATGTAATCGGAGAAGGTGGT"
               "TACGGGATTGTGTATCGTGGCATTTTAACTGATGGAACCAAAGTCGCCGTCAAGAACTTG"
               "CTTAACAATAGGTACTTCCAAATCTTCAATTTTGATTCTAAAGATTGGTCCTTTTACTCT"
               "GTTTCTCAATTTGAGTTTTAGGTATTCTTTGATTTTGTATTGGTTTCATTCTAAATATTC"
               "ATCCTTTACTCAACTTCTAGATAAGGGATTTAGGTATTCTCAAATTTCCGATTTGATTCC"
               "TTTACTCGTTTCTAGATTGGGGTTTTAGGAATTACCAGTTGGGGGTTTTGCAATTTGCGT"
               "AATCAAAGAATTTTATTTGTTGTATTGCTTGGTATTGAAGTTTGTCTCTGTTTCTCTACC"
               "TCGTCATGTAATGTGCTTAGATCCATTAAGTAAATGCTTGTGGATATTTATGTAGATGGT"
               "TAAGAGTGATCGTGATCAGAGTCCTTCTCTTATTTAACTGCATTGCCTGTGAGTTGTGGT"
               "CCTGAAGGTTGTTGTTATTATTGAATTCTATGTATGTATAGATTATGTCATTGGTCTCAT"
               "GTGGTTTTTATGGGTAACGTCTTTACTAATAATAGCACTATGCTTCTGGATTTTGATCTA"
               "TGTGATCTGTAACATTTCTAGTTGGTGTGTCTTTGATTGCCAGGGGTCAAGCAGAGAAGG"
               "AATTCAAAGTAGAAGTGGAAGTCATTGGGCGAGTACGACACAAGAATCTTGTTAGGCTTT"
               "TAGGGTATTGCGTGGAAGGTGCATACAGGTATATATACCGCTGCTCGTATCTCTTTTCCG"
               "GTGTTACAAAAGCGATGTCGTGACCTAATGCTGGGTTCGTTACTATAGGATGCTCGTGTA"
               "TGACTTTGTCGACAATGGTAATTTGGAGCAATGGATTCACGGTGATGTTGGCGATGTCAG"
               "CCCGCTAACTTGGGATATACGTATGAATATTATACTGGGGATGGCCAAAGGGTAAGTCTG"
               "GAATAGCTTTTGAGTTGTCCTCTATGTTTATAAGCTATTGTTGTGTGTAAACCTTTGTTA"
               "TATCTGTAGATTGGCGTATCTACACGAGGGTCTTGAGCCAAAAGTTGTTCATCGGGATAT"
               "TAAATCAAGCAATATCTTACTTGATCGCCAATGGAATGCTAAGGTTTCGGATTTTGGACT"
               "TGCTAAGCTCTTGGGGTCTGAGAGCAGTTATGTGACTACTCGTGTGATGGGAACTTTCGG"
               "GTAAACTATTAAACTCATTAACTCTCTCCTGCAATCTGCAAGGCAGTCTTTAGGAATGTG"
               "AATATTAGGAAATAACTTTTACTTTGTGGGTTGATTTGTTTAGTTATGTAGCACCAGAAT"
               "ACGCTTGCACCGGAATGTTAAACGAGAAGAGTGATATCTATAGCTTCGGAATACTAATCA"
               "TGGAGATAATCACTGGAAGAAACCCGGTTGATTATAGTCGGCCTCAAGGAGAGGTATTTG"
               "ATAAGCATATTCAATCCTCTCTATGTTTTTGTAAATGGTCTTACTATGTGTCGTGGCTGT"
               "AGACAAATCTAGTGGATTGGCTTAAATCAATGGTGGGAAACCGAAGATCAGAAGAAGTTG"
               "TTGATCCGAAAATACCCGAACCACCATCCTCAAAGGCTCTTAAACGAGTGTTGCTTGTAG"
               "CTTTGCGTTGTGTGGATCCTGATGCGAACAAGAGACCTAAAATGGGTCATATCATACATA"
               "TGCTTGAAGCCGAAGATCTACTCTATCGCGATGTATACAAAACAATTTGCCTTTACGTTT"
               "TTACATTTCTTTAAGAGTTTGAAACATGTCTAAAGCTGGGATAATATTTTTGCAGGAACG"
               "CCGAACAACAAGGGACCATGGAAGCCGCGAGAGACAAGAGACAGCTGTGGTTGCTGCCGG"
               "TAGTGAAAGTGGTGAGAGCGGTTCACGGCATCATCAGCAAAAGCAAAGATGAAAAAAGAG"
               "AGTCACTTGGGTTAAGTGATTTCCACACGACCATTATATTTCATTTTTTTTTCCGGTATA"
               "TTATTGTCTCACTTACATAATATTTTCATTCTTTGGCTCTTCGAGTCTTGGCCTTACGAA"
               "ATTAACATTCTTCTTCTTTGTAACCTTTTGTCAACGATCACATTCACATGTCACTTATTA"
               "GATATTGTAATATGTAATGTTTGGACCGACGTCGAAGCATAAGCAGGTGATTGGTCGATG"
               "GATC"),
              ("27488.m000034.TIGRR0.1", 186,
               "TGGTTATGAAAATAATAAGAAATTTAAAGTTAAATGTCTACATGAACATGTTGAAAATCA"
               "CCTTAAGACTGTAAAAAATGAACGAGCAACTATTTCAACACTTCAAGCTAAAAATGGATT"
               "TAGATGGGATGATACTGCAACAAAATATGTTTATGATGAACAAATGATGGTAAGCATTCT"
               "TAATAA")]

    econtroller, lcontroller = seq_len_controllers('pre-mRNA')
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 3
    compare2controllers(expect, econtroller, lcontroller)
    # and now just the longest (1st & 3rd)
    econtroller, lcontroller = seq_len_controllers('pre-mRNA', longest=True)
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 2
    expect = expect[:1] + expect[-1:]
    compare2controllers(expect, econtroller, lcontroller)


def test_get_CDS():
    """test production of CDS sequence, ignoring phase"""
    expect = [("AT1G01540.2.TAIR10", 1419,
               "ATGTCGGTGTACGACGCTGCTTTCCTTAATACAGAGCTTTCGAAACCGACATCGATTTTTGGTCTCCGGC"
               "TATGGGTCGTGATCGGAATCTTACTTGGATCTCTAATTGTCATCGCACTCTTTCTTCTCTCCCTCTGCTT"
               "AACTTCTCGCCGGAAAAATCGAAAGCCGAGAGCCGATTTCGCCTCCGCCGCCATCGCTACACCGCCGATT"
               "TCAAAGGAGATTAAAGAGATCGTTCCGGCGCAGAATCAGTCTGTTCCGGCGGAGATCCAGGTCGATATCG"
               "GGAAGATCGAGCATCGAGTGGTGTTTTCAGATCGAGTGTCGAGTGGTGAGAGTAGAGGAACAGCGAGTGC"
               "AAGTGAAACGGCGTCGTATTCCGGTAGCGGGAATTGTGGGCCGGAGGTGTCGCATCTTGGATGGGGCCGA"
               "TGGTATACTCTGAGAGAGCTTGAAGCGGCCACGAATGGGCTTTGTGAAGAGAATGTAATCGGAGAAGGTG"
               "GTTACGGGATTGTGTATCGTGGCATTTTAACTGATGGAACCAAAGTCGCCGTCAAGAACTTGCTTAACAA"
               "TAGGGGTCAAGCAGAGAAGGAATTCAAAGTAGAAGTGGAAGTCATTGGGCGAGTACGACACAAGAATCTT"
               "GTTAGGCTTTTAGGGTATTGCGTGGAAGGTGCATACAGGATGCTCGTGTATGACTTTGTCGACAATGGTA"
               "ATTTGGAGCAATGGATTCACGGTGATGTTGGCGATGTCAGCCCGCTAACTTGGGATATACGTATGAATAT"
               "TATACTGGGGATGGCCAAAGGATTGGCGTATCTACACGAGGGTCTTGAGCCAAAAGTTGTTCATCGGGAT"
               "ATTAAATCAAGCAATATCTTACTTGATCGCCAATGGAATGCTAAGGTTTCGGATTTTGGACTTGCTAAGC"
               "TCTTGGGGTCTGAGAGCAGTTATGTGACTACTCGTGTGATGGGAACTTTCGGTTATGTAGCACCAGAATA"
               "CGCTTGCACCGGAATGTTAAACGAGAAGAGTGATATCTATAGCTTCGGAATACTAATCATGGAGATAATC"
               "ACTGGAAGAAACCCGGTTGATTATAGTCGGCCTCAAGGAGAGACAAATCTAGTGGATTGGCTTAAATCAA"
               "TGGTGGGAAACCGAAGATCAGAAGAAGTTGTTGATCCGAAAATACCCGAACCACCATCCTCAAAGGCTCT"
               "TAAACGAGTGTTGCTTGTAGCTTTGCGTTGTGTGGATCCTGATGCGAACAAGAGACCTAAAATGGGTCAT"
               "ATCATACATATGCTTGAAGCCGAAGATCTACTCTATCGCGATGAACGCCGAACAACAAGGGACCATGGAA"
               "GCCGCGAGAGACAAGAGACAGCTGTGGTTGCTGCCGGTAGTGAAAGTGGTGAGAGCGGTTCACGGCATCA"
               "TCAGCAAAAGCAAAGATGA"),
              ("AT1G01540.1.TAIR10", 1161,
               "ATGTCGGTGTACGACGCTGCTTTCCTTAATACAGAGCTTTCGAAACCGACATCGATTTTTGGTCTCCGGC"
               "TATGGGTCGTGATCGGAATCTTACTTGGATCTCTAATTGTCATCGCACTCTTTCTTCTCTCCCTCTGCTT"
               "AACTTCTCGCCGGAAAAATCGAAAGCCGAGAGCCGATTTCGCCTCCGCCGCCATCGCTACACCGCCGATT"
               "TCAAAGGAGATTAAAGAGATCGTTCCGGCGCAGAATCAGTCTGTTCCGGCGGAGATCCAGGTCGATATCG"
               "GGAAGATCGAGCATCGAGTGGTGTTTTCAGATCGAGTGTCGAGTGGTGAGAGTAGAGGAACAGCGAGTGC"
               "AAGTGAAACGGCGTCGTATTCCGGTAGCGGGAATTGTGGGCCGGAGGTGTCGCATCTTGGATGGGGCCGA"
               "TGGTATACTCTGAGAGAGCTTGAAGCGGCCACGAATGGGCTTTGTGAAGAGAATGTAATCGGAGAAGGTG"
               "GTTACGGGATTGTGTATCGTGGCATTTTAACTGATGGAACCAAAGTCGCCGTCAAGAACTTGCTTAACAA"
               "TAGGGGTCAAGCAGAGAAGGAATTCAAAGTAGAAGTGGAAGTCATTGGGCGAGTACGACACAAGAATCTT"
               "GTTAGGCTTTTAGGGTATTGCGTGGAAGGTGCATACAGGATGCTCGTGTATGACTTTGTCGACAATGGTA"
               "ATTTGGAGCAATGGATTCACGGTGATGTTGGCGATGTCAGCCCGCTAACTTGGGATATACGTATGAATAT"
               "TATACTGGGGATGGCCAAAGGATTGGCGTATCTACACGAGGGTCTTGAGCCAAAAGTTGTTCATCGGGAT"
               "ATTAAATCAAGCAATATCTTACTTGATCGCCAATGGAATGCTAAGGTTTCGGATTTTGGACTTGCTAAGC"
               "TCTTGGGGTCTGAGAGCAGTTATGTGACTACTCGTGTGATGGGAACTTTCGGTTATGTAGCACCAGAATA"
               "CGCTTGCACCGGAATGTTAAACGAGAAGAGTGATATCTATAGCTTCGGAATACTAATCATGGAGATAATC"
               "ACTGGAAGAAACCCGGTTGATTATAGTCGGCCTCAAGGAGAGGTATTTGATAAGCATATTCAATCCTCTC"
               "TATGTTTTTGTAAATGGTCTTACTATGTGTCGTGGCTGTAG"),
              ("27488.m000034.TIGRR0.1", 186,
               "TGGTTATGAAAATAATAAGAAATTTAAAGTTAAATGTCTACATGAACATGTTGAAAATCACCTTAAGACT"
               "GTAAAAAATGAACGAGCAACTATTTCAACACTTCAAGCTAAAAATGGATTTAGATGGGATGATACTGCAA"
               "CAAAATATGTTTATGATGAACAAATGATGGTAAGCATTCTTAATAA"
               )]
    econtroller, lcontroller = seq_len_controllers('CDS')
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 3
    for item in lcontroller.export_ranges:
        print(item)
    compare2controllers(expect, econtroller, lcontroller)
    # and now just the longest (1st & 3rd)
    econtroller, lcontroller = seq_len_controllers('CDS', longest=True)
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 2
    expect = expect[:1] + expect[-1:]
    compare2controllers(expect, econtroller, lcontroller)


def test_get_spliced_UTRs():
    expect = [('Chr1:195000-199000:780-979', 200,
               "AAGCCTTTCTCTTTAAATTCGTTATCGTTTTTTTTATTTTATCAATTTAATCTTTTTATT"
               "TGTTTTGTTCTTCCTCGATTCAACACTCGATGCTGTGACAAAGCCCAGATTTTGGCCGGA"
               "AAGTTTTGTGTGTTTCCGGCGAAGAATGCGAAGAAACCGGAAGAATCTGTAACGGAATCT"
               "AAGCTAAAAGTTAAAGTACG"),
              ('Chr1:195000-199000:3384-3684', 301,
               "AAAAAGAGAGTCACTTGGGTTAAGTGATTTCCACACGACCATTATATTTCATTTTTTTTT"
               "CCGGTATATTATTGTCTCACTTACATAATATTTTCATTCTTTGGCTCTTCGAGTCTTGGC"
               "CTTACGAAATTAACATTCTTCTTCTTTGTAACCTTTTGTCAACGATCACATTCACATGTC"
               "ACTTATTAGATATTGTAATATGTAATGTTTGGACCGACGTCGAAGCATAAGCAGGTGATT"
               "GGTCGATGGATCAAACTTATCGCTTTGGACGAAATGGACCACAATGATTCTTTTTTAGCT"
               "C"),
              ('Chr1:195000-199000:812-979', 168,
               "TTTATTTTATCAATTTAATCTTTTTATTTGTTTTGTTCTTCCTCGATTCAACACTCGATG"
               "CTGTGACAAAGCCCAGATTTTGGCCGGAAAGTTTTGTGTGTTTCCGGCGAAGAATGCGAA"
               "GAAACCGGAAGAATCTGTAACGGAATCTAAGCTAAAAGTTAAAGTACG"),
              ('', 579,
               "ACAAATCTAGTGGATTGGCTTAAATCAATGGTGGGAAACCGAAGATCAGAAGAAGTTGTT"
               "GATCCGAAAATACCCGAACCACCATCCTCAAAGGCTCTTAAACGAGTGTTGCTTGTAGCT"
               "TTGCGTTGTGTGGATCCTGATGCGAACAAGAGACCTAAAATGGGTCATATCATACATATG"
               "CTTGAAGCCGAAGATCTACTCTATCGCGATGAACGCCGAACAACAAGGGACCATGGAAGC"
               "CGCGAGAGACAAGAGACAGCTGTGGTTGCTGCCGGTAGTGAAAGTGGTGAGAGCGGTTCA"
               "CGGCATCATCAGCAAAAGCAAAGATGAAAAAAGAGAGTCACTTGGGTTAAGTGATTTCCA"
               "CACGACCATTATATTTCATTTTTTTTTCCGGTATATTATTGTCTCACTTACATAATATTT"
               "TCATTCTTTGGCTCTTCGAGTCTTGGCCTTACGAAATTAACATTCTTCTTCTTTGTAACC"
               "TTTTGTCAACGATCACATTCACATGTCACTTATTAGATATTGTAATATGTAATGTTTGGA"
               "CCGACGTCGAAGCATAAGCAGGTGATTGGTCGATGGATC"
               )
    ]
    econtroller, lcontroller = seq_len_controllers('UTR')
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 4
    compare2controllers(expect, econtroller, lcontroller)
    # and now just the longest (1st & 3rd)
    expect = expect[:2]
    econtroller, lcontroller = seq_len_controllers('UTR', longest=True)
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 2
    compare2controllers(expect, econtroller, lcontroller)
    # and -strand UTRs as well (exonexoncds test data)
    expect = [('Gm01b:3312-3460_Gm01b:2196-2196', 150,
               "CTCTGTCTGGACACTTCTAGCATGTAGTAAAAGATATACATATCCAACAAATTGTTAAAT"
               "TTTAAATAATGTGTTATTTGTTTTTTAAAAACTATTTTCAGTTTTCAGCTAAAAAGGGAA"
               "AAAAAAAACTGGCTTTGGATATTTTCCTGG"
               ),
              ('Gm01b:336-479', 144,
               "AGTGATACACTACTACACAACATCATTGTATTTGTACCTTTTTTTTTAGTTGTTTTCTTG"
               "TGTTTGTTTCGGCAAATAGAAGTTGATGTACATGGCATGTCTGATGTTTGTATTTTGTAC"
               "CTAGAATAATGGAAACAATGCGTT"
               ),
               ('Gm01:923-1415_Gm01:2042-2081', 533,
               "TAACGGCGAAAACTTTGTGTCCACCGCCCAGCCCCTCGGCCAATCCCCCAAAACACAAAA"
               "AACTGTTTTTAAAACATAAAAAAAAAACTCATAACATATGAATAATAATAACAATAAAAA"
               "CTAAGAAGCAATAATTATTATTTAATTGTGCATTAAGATATGATTTAAGGGAGATAAGGG"
               "GCTTGTAGAGAGGGTCCTTTCTTTCAGATCTTGCTCCAAAATCAGAACTACCAAAGCTTC"
               "CATCTGATTTGCAGAGAGTTTTTTTTTTTTTTTTTCTTCTCAAGGGACCAACAGAGGCAA"
               "CCAAGGACACTGTAGAGAGAGGAAAGTGGTTCTGTTCTTTCTGTGAGTAACATAACATGC"
               "AATAGATACTATTATGTATAGAGAGAGAACAGAAGAGTTGAAGAAAGATTGTTGTTCAGT"
               "TCGAAAAATTTGCCCATCAAGTACCCTTTTTCTGTTTTCATTTTTTGGTTGGTTTTGCTT"
               "CCTATCCAACCCTGTTGTGGCGTTCCTGGATCGTTCTTCTCTGTGGCGGAGAC"
               ),
              ('Gm01:3696-4056', 361,
               "AGCATTGAATTTTCAAAACACATTTTGATGCTTCGAGTCCAAGGAGTTGCACAACATTGA"
               "TTGATCCATTGACACAGGATTCTTGCTAATTGGTACTTTGCCCCCCTTATTTCTTTCTGA"
               "TATTTTTCTTTCAAAGGTTGGGGATGATGGGGACAAGATTCAGATATATTATTCAAGATT"
               "AGCTGAAAAGTTTTCTGGGGAGGAGCTCTTTTGTCCTTTTTTTTGTTTTTTTTTCTTCCT"
               "TTTATGTTTAAAATTTCAACCACTATTTTGTTACATTTAAATTGGCATCTTCCCCCATTT"
               "CCATT")
              ]
    econtroller, lcontroller = seq_len_controllers('UTR', db=EXONEXONCDS_DB)
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 4
    compare2controllers(expect, econtroller, lcontroller)


def test_get_json_feature():
    controller = JsonExportController(db_path_in='sqlite:///' + EXPORTING_DB)
    f = controller.session.query(orm.Feature).filter(orm.Feature.type == types.GEENUFF_CDS).first()
    print(f.type)
    t = f.transcript_pieces[0].transcript
    coord = controller.session.query(orm.Coordinate).first()
    fh = FeatureJsonable(data=f)
    th = TranscriptJsonable(data=t)
    print(fh.to_jsonable(f, coord, 790, 3000, True, transcript=t))
    #print(json.dumps(th.to_jsonable(t, coord, 790, 3000, True), indent=2))
    print('------')
    sls = controller.session.query(orm.SuperLocus).all()
    #print(len(sls))
    for sl in sls:
        slh = SuperLocusJsonable(sl)
        #print(json.dumps(slh.to_jsonable(sl, coord, 790, 3000, True), indent=2))
    #print(controller.session.query(orm.Genome).all())
    meh = controller.coordinate_range_to_jsonable('dummy', seqid='Chr1:195000-199000', start=1, end=3900, is_plus_strand=True)
    assert len(meh) == 1
    assert meh[0]['coordinate_piece']['seqid'] == 'Chr1:195000-199000'
    assert len(meh[0]['coordinate_piece']['sequence']) == 3899
    assert len(meh[0]['super_loci']) == 1
    print(json.dumps(meh, indent=2))
    # todo, slightly more thorough testing


def helper_get_species(pks, session):
    out = []
    for pk in pks:
        seqid = session.query(orm.Genome).filter(orm.Genome.id == pk).first().species
        out.append(seqid)
    return out


def test_exporter_all_or_1_transcript():
    controller = GeenuffExportController(db_path_in='sqlite:///' + EXPORTING_DB)
    coords = controller.genome_query(all_transcripts=True)
    assert len(coords.keys()) == 2
    # first coordinate has 1 super locus with two transcripts. pk=1, len=4000
    features_c1 = coords[(1, 4000)]
    assert len([x for x in features_c1 if x.type.value == types.GEENUFF_TRANSCRIPT]) == 2
    # second coord has 1 sl with 1 transcript
    features_c2 = coords[(2, 3000)]
    assert len([x for x in features_c2 if x.type.value == types.GEENUFF_TRANSCRIPT]) == 1
    # queried with longest both should have only one transcript remaining
    coords = controller.genome_query(all_transcripts=False)
    assert len(coords.keys()) == 2
    for coord, features in coords.items():
        assert len([x for x in features if x.type.value == types.GEENUFF_TRANSCRIPT]) == 1


def test_exporter_return_super_loci():
    # this test is here to make sure I don't break something I didn't notice while refactoring
    # in the long run, I almost certainly want to _change_ the behaviour so genome_query does
    # not have two different output types.
    controller = GeenuffExportController(db_path_in='sqlite:///' + EXPORTING_DB)
    genome_coords = controller.genome_query(all_transcripts=True, return_super_loci=True)
    sl = genome_coords[0][0]
    assert isinstance(sl, orm.SuperLocus)
    coord_id = genome_coords[1][1]
    assert coord_id == "27488:1-3000"
