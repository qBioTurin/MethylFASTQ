#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, argparse
import os.path
import string, re
from Bio import SeqIO
from multiseq import ChromoSeq
from timeit import default_timer as timer

from seqparams import SeqParams

# TODO
# - percentuale metilazione - ok
# - dove mappano le read (su quali citosine)
# - confronto allineamento con dataset simulato \


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MethylSeq")

    parser.add_argument("--in", dest="fasta_file", metavar="fasta-file",
                        action="store", type=str, required=True,
                        help="Path of FASTA file containing the genome to sequencing.")
    parser.add_argument("--out", dest="output_path", metavar="output-file-path",
                        action="store", type=str, required=True,
                        help="Path of output files (.fastq and .cpg)")
    parser.add_argument("--seq", dest="seq_mode",
                        choices=["single_end", "paired_end"], default="single_end",
                        help="Sequencing type")
    parser.add_argument("--lib", dest="lib_mode",
                        choices=["directional", "non_directional"], default="directional",
                        help="")
    parser.add_argument("--chr", dest="chr", metavar="chromosome-id",
                        nargs="*", action="store", type=str,
                        help="Chromosomes to sequencing")
    parser.add_argument("--fragment", dest="fsize", metavar="fragment-size",
                        action="store", type=int, default=200,
                        help="Size of fragments in the fragmenation step")
    parser.add_argument("--read", dest="rlength", metavar="read-length",
                        action="store", type=int, default=200,
                        help="")
    parser.add_argument("--processes", dest="processes", metavar="num-processes",
                        action="store", type=int, default=2,
                        help="")
    parser.add_argument("--temp", dest="temp_dir", metavar="temporary-directory",
                        action="store", type=str, help="")
    parser.add_argument("--cg", dest="p_cg", metavar="CG-methylation-probability",
                        action="store", type=float, default=0.82,
                        help="")
    parser.add_argument("--chg", dest="p_chg", metavar="CHG-methylation-probability",
                        action="store", type=float, default=0.01,
                        help="")
    parser.add_argument("--chh", dest="p_chh", metavar="CHH-methylation-probability",
                        action="store", type=float, default=0.02,
                        help="")
    parser.add_argument("--coverage", dest="coverage", metavar="coverage",
                        action="store", type=float, default=6.3,
                        help="")
    parser.add_argument("--snp", dest="snp", metavar="snp-probability",
                        action="store", type=float, default=0,
                        help="")
    parser.add_argument("--error", dest="error_rate", metavar="sequencing-error-probability",
                        action="store", type=float, default=0,
                        help="")
    # parser.add_argument("--config", dest="config", metavar="configuration-file",
    #                     action="store", type=str, default=False, #?
    #                     help="")

    args = parser.parse_args()
    params = None
    try:
        params = SeqParams().load_params(args)
    except Exception as e:
        print(e)
        sys.exit("Aborted.")

    params.params_summary()

    regex = re.compile('[%s]' % re.escape(string.punctuation))
    meth_probs = (params["p_meth_cg"], params["p_meth_chg"], params["p_meth_chh"])

    for fasta_record in SeqIO.parse(params["input_file"], "fasta"):
        t, flag, chromo = fasta_record.description.split(), True, None

        try:
            chromo = regex.sub("", t[t.index("chromosome") + 1])
        except ValueError:
            flag = False

        if params["chromosomes"] is None or (flag and chromo in params["chromosomes"]):
            print("Sequenzio il cromosoma {}: {} ({} bp)".format(chromo, fasta_record.id, len(fasta_record)), flush=True)
            print(fasta_record.description)

            ChromoSeq(fasta_record, meth_probs).sequencing(params)
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
