#! /usr/bin/env python3
import argparse

from geenuff.applications.exporters.sequence import FastaExportController

def main(args):
    controller = FastaExportController(args.db_path_in, args.longest)
    if args.mode == "introns":
        controller.prep_intron_exports(args.genomes, args.exclude_genomes)
    elif args.mode == "pre-mRNA":
        controller.prep_transcript_exports(args.genomes, args.exclude_genomes)
    else:
        raise NotImplementedError("I lied, only mode=introns, pre-mRNA is implemented so far, not {}".format(args.mode))
    coords_ids = controller._get_coords_by_genome_query(args.genomes, args.exclude_genomes)
    #x = controller.get_super_loci_by_coords(coords_ids)
    controller.write_fa(args.out)


if __name__ == "__main__":
    implemented_modes = ["mRNA", "pre-mRNA", "CDS", "introns", "exons", "UTR"]
    parser = argparse.ArgumentParser()
    parser.add_argument('--db-path-in', type=str, required=True,
                        help='Path to the Geenuff SQLite input database.')
    parser.add_argument('-m', '--mode', type=str, required=True,
                        help="which type of sequence to export, one of {}".format(implemented_modes))
    parser.add_argument('-o', '--out', type=str,
                        help="output fasta file path, (default is stdout)")
    parser.add_argument('--genomes', type=str, default='',
                        help='Comma separated list of species names to be exported. '
                             'If empty all genomes in the db are used.')
    parser.add_argument('--exclude-genomes', type=str, default='',
                        help='Comma separated list of species names to be excluded. '
                             'If empty all genomes in the db are used.')

    parser.add_argument('-l', '--longest', action="store_true",
                        help="ignore all but the longest transcript per gene")
    args = parser.parse_args()

    main(args)
