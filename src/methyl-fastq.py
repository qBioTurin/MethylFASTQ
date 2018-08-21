#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# pnr: g788g5 7-12a


import argparse, sys
import os, pathlib
import random, string
import csv, re
from Bio import SeqIO
from multiseq import ChromoSeq
from timeit import default_timer as timer

#from seqparams import SeqParams

# TODO
# - percentuale metilazione - ok
# - dove mappano le read (su quali citosine)
# - confronto allineamento con dataset simulato \

class MethylFASTQ(object):
    def __init__(self, args):
        self.params = self.__parse_args(args)
        self.regions = self.__load_regions(args.regions)


    def run(self):
        print("Chosen mode: {}".format("targeted" if self.regions else "WGBS"))

        self.__run_targeted() if self.regions else self.__run_wgbs()

        #rimuovo directory dei file temporanei
        print("Rimuovo directory {}".format(args.temp_dir))
        os.rmdir(args.temp_dir)

    def __parse_args(self, args):
        if args.temp_dir is None:
            temp_dir = os.fspath("{}/.cache/methylfastq/".format(pathlib.Path.home()))

            #creo directory principale, se non presente
            if not os.path.exists(temp_dir):
                print("MethylFASTQ temporary files will be written in {}".format(temp_dir))
                os.makedirs(temp_dir)

            #creo directory per la run corrente
            while True:
                random_stuff = "".join(random.choices(string.ascii_lowercase, k=2))
                if not os.path.exists(temp_dir + random_stuff):
                    temp_dir += random_stuff + "/"
                    os.makedirs(temp_dir)
                    print("Temporary files of current run will be written in {}".format(temp_dir))
                    break

            args.temp_dir = temp_dir

        if args.fragment_size == 0:
            args.fragment_size = args.read_length

        return args

    def __run_wgbs(self):
        regex = re.compile('[%s]' % re.escape(string.punctuation))

        for fasta_record in SeqIO.parse(self.params.fasta_file, "fasta"):
            descr = fasta_record.description.split()
            flag, chromo = True, None

            try:
                chromo = regex.sub("", descr[descr.index("chromosome") + 1])
            except ValueError:
                flag = False

            if self.params.chr is None or (flag and chromo in self.params.chr):
                print("Sequencing {}: {} ({} bp)\n{}".format(chromo, fasta_record.id, len(fasta_record), fasta_record.description), flush=True)

                ChromoSeq(fasta_record).sequencing(self.params)

    def __run_targeted(self):
        for fasta_record in SeqIO.parse(self.params.fasta_file, "fasta"):
            if fasta_record.id in self.regions:
                curr_chr = self.regions[fasta_record.id]
                print("{} --> {}".format(fasta_record.id, curr_chr))

                ChromoSeq(fasta_record, target_regions=curr_chr).sequencing(self.params)

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

    #I/O parameters: input, output and temporary directories
    parser.add_argument("-i", "--in", dest="fasta_file", metavar="fasta-file",
                        action="store", type=str, required=True,
                        help="Path of FASTA file containing the genome to sequencing.")
    parser.add_argument("-o", "--out", dest="output_path", metavar="output-file-path",
                        action="store", type=str, required=True,
                        help="Path of output files (.fastq and .cpg)")
    parser.add_argument("--temp", dest="temp_dir", metavar="temporary-directory",
                        action="store", type=str, help="")
    #sequencing mode and library mode
    parser.add_argument("--seq", dest="seq_mode",
                        choices=["single_end", "paired_end"], default="single_end",
                        help="Sequencing type")
    parser.add_argument("--lib", dest="lib_mode",
                        choices=["directional", "non_directional"], default="directional",
                        help="")
    #chromosomes to be sequenced
    parser.add_argument("--chr", dest="chr", metavar="chromosome-id",
                        nargs="*", action="store", type=str,
                        help="Chromosomes to be sequence")
    parser.add_argument("--regions", dest="regions", metavar="target-regions",
                        action="store", type=str,
                        help="Target regions to ")
    #coverage
    parser.add_argument("--coverage", dest="coverage", metavar="coverage",
                        action="store", type=float, default=6.3,
                        help="")

    #fragment and read length
    parser.add_argument("--fragment", dest="fragment_size", metavar="fragment-size",
                        action="store", type=int, default=0, #default = read_length
                        help="Size of fragments in the fragmenation step")
    parser.add_argument("--read", dest="read_length", metavar="read-length",
                        action="store", type=int, default=150,
                        help="")
    #parallel option: num of processes to be used
    parser.add_argument("--processes", dest="num_processes", metavar="num-processes",
                        action="store", type=int, default=2,
                        help="")
    #methylation probabilities
    parser.add_argument("--cg", dest="p_cg", metavar="CG-methylation-probability",
                        action="store", type=float, default=1.0,
                        help="")
    parser.add_argument("--chg", dest="p_chg", metavar="CHG-methylation-probability",
                        action="store", type=float, default=1.0,
                        help="")
    parser.add_argument("--chh", dest="p_chh", metavar="CHH-methylation-probability",
                        action="store", type=float, default=1.0,
                        help="")
    #sequencing error and snp probabilities
    parser.add_argument("--snp", dest="snp_rate", metavar="snp-probability",
                        action="store", type=float, default=0,
                        help="")
    parser.add_argument("--error", dest="error_rate", metavar="sequencing-error-probability",
                        action="store", type=float, default=0,
                        help="")

    args = parser.parse_args()

    print(args.fragment_size)

    MethylFASTQ(args).run()


#
#
#     if False and args.regions:
#         target_regions = dict()
#
#         #load .bed file
#         try:
#             with open(args.regions) as f:
#                 for line_number, line in enumerate(csv.reader(f, delimiter="\t")):
#                     try:
#                         chromosome = line[0]
#                         region = (int(line[1]), int(line[2]))
#
#                         if chromosome not in target_regions:
#                             target_regions[chromosome] = list()
#
#                         target_regions[chromosome].append(region)
#                     except IndexError:
#                         print("Malformed line num {} --> {}. Skipped".format(line_number+1, line))
#                     except ValueError:
#                         print("Error in line num {}: second and third fields must be integers. Skipped".format(line_number+1, line))
#
# #            print(target_regions)
#
#             for fasta_record in SeqIO.parse(args.fasta_file, "fasta"):
#                 if fasta_record.id in target_regions:
#                     print("{} --> {}".format(fasta_record.id, target_regions[fasta_record.id]))
#                     regions = target_regions[fasta_record.id]
#                     ChromoSeq(fasta_record, target_regions=regions).sequencing(args)
# #                    print(fasta_record.id)
#                     #fai qualcosa
#
#         except FileNotFoundError:
#             sys.exit("File {} does not exist. Aborted.".format(args.regions))
#
#     else:
#         regex = re.compile('[%s]' % re.escape(string.punctuation))
#
#         for fasta_record in SeqIO.parse(args.fasta_file, "fasta"):
# #            print(dir(fasta_record))
# #            print(fasta_record.id)
#             t, flag, chromo = fasta_record.description.split(), True, None
#
#             try:
#                 chromo = regex.sub("", t[t.index("chromosome") + 1])
#             except ValueError:
#                 flag = False
#
#             if args.chr is None or (flag and chromo in args.chr):
#                 print("Sequenziooooooooooooooooo il cromosoma {}: {} ({} bp)".format(chromo, fasta_record.id, len(fasta_record)), flush=True)
#                 print(fasta_record.description)
#
#                 ChromoSeq(fasta_record).sequencing(args)




#    sys.exit("Grazie e arrivederci")





    #
    #
    #
    # regex = re.compile('[%s]' % re.escape(string.punctuation))
    # meth_probs = (params["p_meth_cg"], params["p_meth_chg"], params["p_meth_chh"])
    #
    # for fasta_record in SeqIO.parse(params["input_file"], "fasta"):
    #     t, flag, chromo = fasta_record.description.split(), True, None
    #
    #     try:
    #         chromo = regex.sub("", t[t.index("chromosome") + 1])
    #     except ValueError:
    #         flag = False
    #
    #     if params["chromosomes"] is None or (flag and chromo in params["chromosomes"]):
    #         print("Sequenzio il cromosoma {}: {} ({} bp)".format(chromo, fasta_record.id, len(fasta_record)), flush=True)
    #         print(fasta_record.description)
    #
    #         ChromoSeq(fasta_record).sequencing(params)
    #
    #        multi = ChromoSeq(fasta_record, meth_probs=(params.p_cg, params.p_chg, params.p_chh))
    #        multi.sequencing(params.processes, params.seq_mode, params.lib_mode, params.fsize, params.rlength, params.output_path, params.temp_dir)


    ##############################################################################
    #         start = timer()
    #         chr_seq = ChromosomeSequencing(params.output_path, fasta_record)
    # #        chr_seq.load_methylation(params.methyl) #se methyl è None non ha effetto
    #         chr_seq.sequencing(params.seq_mode, params.lib_mode, params.fsize, params.rlength)
    #
    #         print("Sequencing of {} done in {}".format(fasta_record.id, format_time(timer()-start)), flush=True)
    ##############################################################################
