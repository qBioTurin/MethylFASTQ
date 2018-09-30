#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, sys
import os, pathlib
import random, string
import csv
from Bio import SeqIO
import sequencing as seq

from timeit import default_timer as timer


class MethylFASTQ(object):
    def __init__(self, args):
        self.params = self.__parse_args(args)
        self.regions = self.__load_regions(args.regions)
        self.__stats = seq.Stats()



    def run(self):
        print("Chosen mode: {}".format("targeted" if self.regions else "WGBS"))

        start = timer()

        self.__run_targeted() if self.regions else self.__run_wgbs()

        elapsed = timer() - start

        print("Num reads: {}\nNum C: {}\nNum bases: {}\n".format(self.__stats.nreads, self.__stats.ncytosines, self.__stats.nbases))

        print("Elapsed time: {}".format(seq.format_time(elapsed)))


    def __parse_args(self, args):
        #check esistenza directory output
        if not os.path.exists(args.output_path):
            print("Output directory {} does not exist.".format(args.output_path), end=" ")
            os.makedirs(args.output_path)
            print("Now it does.")

        if args.fragment_size == 0:
            args.fragment_size = args.read_length

        if args.num_processes < 1:
            raise ValueError("Number of processes must be greater than zero.")

        return args

    def __run_wgbs(self):
        num_reads, num_c = 0, 0

        for fasta_record in SeqIO.parse(self.params.fasta_file, "fasta"):
            if self.params.chr is None or fasta_record.id in self.params.chr:
                print("Sequencing {}: {} bp".format(fasta_record.id, len(fasta_record)), flush=True)

                stats = seq.ChromosomeSequencer(fasta_record).sequencing(self.params)
                self.__stats.update(stats)



    def __run_targeted(self):
        for fasta_record in SeqIO.parse(self.params.fasta_file, "fasta"):
            if fasta_record.id in self.regions:
                curr_chr = self.regions[fasta_record.id]
                print("{} --> {}".format(fasta_record.id, curr_chr))

                stats = seq.ChromosomeSequencer(fasta_record, target_regions=curr_chr).sequencing(self.params)
                self.__stats.update(stats)

    def __load_regions(self, filename):
        target_regions = None

        try:
            with open(filename) as f:
                target_regions = dict()
                csvr = csv.reader(f, delimiter='\t')

                for ln, line in enumerate(csvr, 1):
                    try:
                        chromosome = line[0]
                        interval = int(line[1]), int(line[2])

                        if chromosome not in target_regions:
                            target_regions[chromosome] = list()
                        target_regions[chromosome].append(interval)
                    except IndexError:
                        print("Malformed line #{}: {}.\tSkipped".format(ln, line))
                    except ValueError:
                        print("Error in line #{}: {}.\tFields #2 and #3 must have integer values. Skipped".format(ln, line))
        except FileNotFoundError:
            sys.exit("File {} does not exist. Aborted".format(filename))
        except TypeError:
            pass #filename == None -> return None

        return target_regions




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="methyl-fastq")

    #I/O parameters: input, output
    parser.add_argument("-i", "--in", dest="fasta_file", metavar="fasta-file",
                        action="store", type=str, required=True,
                        help="Path of the FASTA file containing the genome to be sequenced.")
    parser.add_argument("-o", "--out", dest="output_path", metavar="output-file-path",
                        action="store", type=str, required=True,
                        help="Path of output files (.fastq & .cpg)")

    #sequencing mode and library mode
    parser.add_argument("--seq", dest="seq_mode",
                        choices=["single_end", "paired_end"], default="single_end",
                        help="Type of reads to produce")

    parser.add_argument("--lib", dest="lib_mode",
                        choices=["directional", "non_directional"], default="directional",
                        help="Library type to produce")
    #chromosomes to be sequenced
    parser.add_argument("--chr", dest="chr", metavar="chromosome-id",
                        nargs="*", action="store", type=str,
                        help="List of chromosomes to be sequenced (WGBS)")

    parser.add_argument("--regions", dest="regions", metavar="target-regions",
                        action="store", type=str,
                        help="Genomic regions to be sequenced (targeted bisulfite sequencing)")
    #coverage
    parser.add_argument("--coverage", dest="coverage", metavar="coverage",
                        action="store", type=float, default=6.3,
                        help="Depth of coverage")

    #fragment and read length
    parser.add_argument("--fragment", dest="fragment_size", metavar="fragment-size",
                        action="store", type=int, default=0, #default = read_length
                        help="Size of fragments in the fragmenation step")

    parser.add_argument("--read", dest="read_length", metavar="read-length",
                        action="store", type=int, default=150,
                        help="Read length of sequencing step")
    #parallel option
    parser.add_argument("--processes", dest="num_processes", metavar="num-processes",
                        action="store", type=int, default=2,
                        help="Number of processes to be used during sequencing step")

    parser.add_argument("--buffer", dest="buffer_size", metavar="buffer-size",
                        action="store", type=int, default=10**6,
                        help="Buffer size of each process during sequencing step")
    #methylation probabilities
    parser.add_argument("--cg", dest="p_cg", metavar="CG-methylation-probability",
                        action="store", type=float, default=1.0,
                        help="Methylation probability in CG context")

    parser.add_argument("--chg", dest="p_chg", metavar="CHG-methylation-probability",
                        action="store", type=float, default=1.0,
                        help="Methylation probability in CHG context")

    parser.add_argument("--chh", dest="p_chh", metavar="CHH-methylation-probability",
                        action="store", type=float, default=1.0,
                        help="Methylation probability in CHH context")

    #sequencing error and snp probabilities
    parser.add_argument("--snp", dest="snp_rate", metavar="snp-probability",
                        action="store", type=float, default=0,
                        help="Probability to set a SNP")

    parser.add_argument("--error", dest="error_rate", metavar="sequencing-error-probability",
                        action="store", type=float, default=0,
                        help="Probability to set a sequencing error")
    #quality of reads
    parser.add_argument("--maxq", dest="max_quality", metavar="max-phred-score",
                        action="store", type=int, default=40,
                        help="Max phred score in the reads")

    parser.add_argument("--minq", dest="min_quality", metavar="min-phred-score",
                        action="store", type=int, default=10,
                        help="Min phred score in the reads (not implemented)") #not implemented yet!!


    args = parser.parse_args()

    MethylFASTQ(args).run()
