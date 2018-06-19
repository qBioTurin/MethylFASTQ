#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
from io import StringIO
import random, csv
import os, pickle, multiprocessing as mp
from timeit import default_timer as timer
import enum

import dna


#TODO:
#1. il nome dei record fastq nel sequenziamento paired_end deve finire con /1 o /2
#2.



def format_time(time):
    strtime = "{:03.5f} m".format(time/60) if time >= 60 else "{:03.5f} s".format(time)
    return strtime

class CytosineContext(enum.Enum):
    CG = 1
    CHG = 2
    CHH = 3

    def __str__(self):
        return ""


class ChromoSeq(object):
    """ Un oggetto ChromoSeq legge un cromosoma da un file FASTA e
    scrive su un file FASTQ le read prodotte dal sequenziamento.
    È possibile scegliere il tipo di library (single-end o paired-end),
    se la library è direzionale o non direzionale, la lunghezza dei frammenti
    e delle read """

    def __init__(self, chromosome, meth_probs):
        """Parserizza il cromosoma individuando gli intervalli da sequenziare e scartando
        i nucleotidi indefiniti (basi N)"""

        self.__chromoId = chromosome.id
        self.__genome_size = len(chromosome.seq)
        self.__fragments = list()

        print("Parsing {} sequence... ".format(chromosome.id), end="", flush=True)
        chromosome = str(chromosome.seq).lower()

        #get salient fragments
        start, i = timer(), 0

        while i < self.__genome_size:
            while i < self.__genome_size and chromosome[i] == 'n':
                i += 1
            begin = i
            while i < self.__genome_size and chromosome[i] != 'n':
                i += 1

            t = (chromosome[begin:i], begin, i)
            self.__fragments.append(t)
        else:
            last = self.__fragments.pop()
            if last[1] < last[2]:
                self.__fragments.append(last)
        print("{} fragments found in {}.".format(len(self.__fragments), format_time(timer() - start)), flush=True)


    def sequencing(self, params):
        """Sequenzia i frammenti individuati dal costruttore parallelizzando il
         lavoro con n processi. """

        self.__fragments = [{"seq": x[0], "from": x[1], "to": x[2], "params": params} for x in self.__fragments]
        seq_function, merge_function = (self.single_end_reads, self.merge_se_file)\
                                        if params.seq_mode == "single_end"\
                                        else (self.paired_end_reads, self.merge_pe_file)

        start = timer()
        with mp.Pool(params.num_processes) as p:
            p.map(seq_function, self.__fragments)
        print("Sequencing of {} completed in {}".format(self.__chromoId, format_time(timer()-start)))

        merge_function(params)
        self.merge_ch3(params)
        print()


    def __get_output_filename(self, params):
        fasta = "".join(params.input_file.split("/")[-1].split(".")[:-1])
        se_pe = "se" if params.seq_mode == "single_end" else "pe"
        dir_nondir = "dir" if params.lib_mode == "directional" else "undir"
        return "{}/{}_{}_f{}r{}_{}".format(params.output_dir, fasta, se_pe, params.fragment_size, params.read_length, dir_nondir)

    def merge_se_file(self, params):
        """ """

        output_filename = self.__get_output_filename(params)
        print("Writing final fastq file...", flush=True, end="")

        with open("{}.fastq".format(output_filename), "a") as out:
            for e in self.__fragments:
                input_name = "{}/{}.seq".format(params.temp_dir, e["from"])
        #        print("-reading {}".format(input_name), flush=True)

                try:
                    with open(input_name, "rb") as input_file:
                        try:
                            while True:
                                record = pickle.load(input_file)
                                SeqIO.write(record, out, "fastq-illumina")
                        except EOFError:
                            print(".", flush=True, end="")
                    os.remove(input_name)
                except FileNotFoundError as e:
                    pass

    def merge_pe_file(self, params):
        """Genera i file FASTQ del sequenziamento paired-end leggendo i file binari temporanei
        prodotti nello step di sequenziamento. """

        output_filename = self.__get_output_filename(params)
        print("Writing final fastq files...", flush=True, end="")

        with open("{}_R1.fastq".format(output_filename), "a") as r1, open("{}_R2.fastq".format(output_filename), "a") as r2:
            for e in self.__fragments:
                input_name = "{}/{}.seq".format(params.temp_dir, e["from"])
#                print("Reading {}...".format(input_name), flush=True)

                try:
                    with open(input_name, "rb") as input_file:
                        try:
                            while True:
                                record = pickle.load(input_file)
                                SeqIO.write(record[0], r1, "fastq-illumina")
                                SeqIO.write(record[1], r2, "fastq-illumina")
                        except EOFError:
                            print(".", flush=True, end="") #done!!
                    os.remove(input_name)
                except FileNotFoundError as e:
                    pass

    def merge_ch3(self, params):
        """Genera un unico file contenente i dati di metilazione delle citosine
        unendo i file relativi ai diversi frammenti prodotti dai vari processi
        nello step di sequenziamento"""

#        output_file = "{}.ch3".format(output_prefix)
        output_filename = self.__get_output_filename(params)
        print("\nWriting final methylation file...", flush=True, end="")

        with open("{}.ch3".format(output_filename), "a") as out:
            csvwriter = csv.writer(out, delimiter=" ")

            for e in self.__fragments:
                input_name = "{}/{}.ch3".format(params.temp_dir, e["from"])

                with open(input_name) as infile:
                    csvreader = csv.reader(infile, delimiter=" ")
                    for row in csvreader:
                        csvwriter.writerow(row)

                os.remove(input_name)
                print(".", flush=True, end="") #done!!


    def single_end_reads(self, data):
        params = data["params"]

        seq = FragmentSeq(self.__chromoId, data["seq"], data["from"], data["to"], self.__genome_size, data["params"])
        seq.single_end_sequencing()
        print(".", flush=True, end="")

    def paired_end_reads(self, data):
        params = data["params"]

        seq = FragmentSeq(self.__chromoId, data["seq"], data["from"], data["to"], self.__genome_size, data["params"])
        seq.paired_end_sequencing()
        print(".", flush=True, end="")


class FragmentSeq(object):
    """..."""

    def __init__(self, chr_id, sequence, begin_seq, end_seq, gsize, seqparams):
        #informazioni sulla sequenza
        self.__chromoId = chr_id
        self.__sequence = sequence
        self.__genome_size = gsize
        #posizione del frammento nel genoma
        self.__offset_forward = begin_seq
        self.__offset_reverse = gsize - end_seq
        #dimensione buffer + parametri vari
        self.__buffer_size = 10**4
        self.__seqparams = seqparams
        self.__p_meth = seqparams.p_meth
        #dati sulle citosine
        self.__cytosines_forward = dict()
        self.__cytosines_reverse = dict()
        ##############################
        self.__set_snp()
        self.__initialize_cytosines()


    ###################### Sequencing ######################

    def fragmentation(self):
        """ """
        max_step = 2*int(round(self.__seqparams.fragment_size / self.__seqparams.coverage))
        i = 0

        while (i + self.__seqparams.fragment_size) < len(self.__sequence):
            fragment = self.__sequence[i: i + self.__seqparams.fragment_size]

            if len(fragment) == self.__seqparams.fragment_size:
                yield dna.fragment(fragment, i, i + self.__seqparams.fragment_size)\
                        .methylate(self.methylate_cytosine)

            i += random.randint(1, max_step) #in media ogni base è coperta da C reads


    def single_end_sequencing(self):
        """Produce le read single-end del frammento e le salva su un file fastq temporaneo"""

        fragment_size = self.__seqparams.fragment_size
        read_length = self.__seqparams.read_length
        directional = self.__seqparams.lib_mode == "directional"
        seq_length = len(self.__sequence)
        error_rate = self.__seqparams.error_rate

        reads = list()

        for fragment in self.fragmentation():
            #metilazione?
            bsfs = fragment.forward_strand.bisulfite()
            fastq_fw = bsfs.single_end_sequencing(read_length)\
                        .set_sequencing_errors(error_rate)\
                        .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

            bsrs = fragment.reverse_strand.bisulfite()
            fastq_rv = bsrs.single_end_sequencing(read_length)\
                        .set_sequencing_errors(error_rate)\
                        .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

            reads.extend([fastq_fw, fastq_rv])

            if not directional:
                fastq_fwrc = bsfs.reverse_complement()\
                            .single_end_sequencing(read_length)\
                            .set_sequencing_errors(error_rate)\
                            .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

                fastq_rvrc = bsrs.reverse_complement()\
                            .single_end_sequencing(read_length)\
                            .set_sequencing_errors(error_rate)\
                            .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

                reads.extend([fastq_fwrc, fastq_rvrc])

            if len(reads) > self.__buffer_size:
                self.__persist_record(reads)
                reads.clear()
        else:
            if len(reads) > 0:
                self.__persist_record(reads)
                reads.clear()

        self.__persist_methylation()


    def paired_end_sequencing(self):
        """Produce le read paired-end del frammento e le salva su un file temporaneo"""

        fragment_size = self.__seqparams.fragment_size
        read_length = self.__seqparams.read_length
        directional = self.__seqparams.lib_mode == "directional"
        seq_length = len(self.__sequence)
        error_rate = self.__seqparams.error_rate

        reads = list()

#        for (fragment, begin, end) in self.__fragmentation():
        for fragment in self.fragmentation():
            bsfs = fragment.forward_strand.bisulfite()
            fastq_fw = bsfs.paired_end_sequencing(read_length)\
                        .set_sequencing_errors(error_rate)\
                        .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

            bsrs = fragment.reverse_strand.bisulfite()
            fastq_rv = bsrs.paired_end_sequencing(read_length)\
                        .set_sequencing_errors(error_rate)\
                        .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

            reads.extend([fastq_fw, fastq_rv])

            if not directional:
                fastq_fwrc = bsfs.reverse_complement()\
                            .paired_end_sequencing(read_length)\
                            .set_sequencing_errors(error_rate)\
                            .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

                fastq_rvrc = bsrs.reverse_complement()\
                            .paired_end_sequencing(read_length)\
                            .set_sequencing_errors(error_rate)\
                            .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

                reads.extend([fastq_fwrc, fastq_rvrc])

            if len(reads) > self.__buffer_size:
                self.__persist_record(reads)
                reads.clear()
        else:
            if len(reads) > 0:
                self.__persist_record(reads)
                reads.clear()

        self.__persist_methylation()

    ###################### Methylation ######################

    def __initialize_cytosines(self):
        """Parserizza il genoma e indicizza le citosine sui due strand"""

#        print("Init fragment from {} to {}".format(self.__offset_forward, self.__genome_size - self.__offset_reverse), flush=True)
        limit = len(self.__sequence)

        for pos, base in enumerate(self.__sequence):
            #strand +
            if base == 'c':
                context = CytosineContext.CHH
                if pos+1 < limit and self.__sequence[pos+1] == 'g':
                    context = CytosineContext.CG
                elif pos+2 < limit and self.__sequence[pos+2] == 'g':
                    context = CytosineContext.CHG
                #initialize new cythosine
                self.__cytosines_forward[pos] = {"context": context, "nmeth": 0, "ntot": 0}
            #strand -
            elif base == 'g':
                context = CytosineContext.CHH
                if pos-1 >= 0 and self.__sequence[pos-1] == "c":
                    context = CytosineContext.CG
                elif pos-2 >= 0 and self.__sequence[pos-2] == "c":
                    context = CytosineContext.CHG
                #initialize new cythosine
                self.__cytosines_reverse[pos] = {"context": context, "nmeth": 0, "ntot": 0}


    def methylate_cytosine(self, base, position):
        strand_dict = self.__cytosines_forward if base == 'c' or base == 'C' else self.__cytosines_reverse
        state = base.upper()
        if position in strand_dict:
            cytosine = strand_dict[position]
            cytosine["ntot"] += 1

            if random.uniform(0, 1) <= self.__p_meth[cytosine["context"]]:
                state = state.lower()
                cytosine["nmeth"] += 1
        return state


    ###################### Introduzione mutazioni ######################

    def __set_snp(self):
        """Setta SNP random sulla reference"""

        self.__sequence = "".join([random.sample("actg", 1)[0] if random.uniform(0, 1) < self.__seqparams.snp_rate else base\
                                   for base in self.__sequence])

    ###################### Record fastq ######################

    def __set_quality(self, sequence):
        """Genera una quality verosimile"""
        quality = [42] * len(sequence)

        return quality

    ###################### Persistenza file temporanei ######################

    def __persist_record(self, record_list):
        """Accoda i record contenenti nella lista @record_list nel file binario @file_name"""

        filename = "{}/{}.seq".format(self.__seqparams.temp_dir, self.__offset_forward)

        with open(filename, "ab") as output:
            for record in record_list:
                pickle.dump(record, output, pickle.HIGHEST_PROTOCOL)

    def __persist_methylation(self):
        """Write a file with a line for each cythosine. Each line contains:
        - chromosome id, position, strand, cythosine context
        - times C is methylated, total times C appears in read file, ratio"""

        ch3_file = "{}/{}.ch3".format(self.__seqparams.temp_dir, self.__offset_forward)

        with open(ch3_file, "w") as of:
            csvfile = csv.writer(of, delimiter="\t")
            prep = [(self.__cytosines_forward, "+", self.__offset_forward),
                    (self.__cytosines_reverse, "-", self.__offset_forward)]
#                    (self.__cytosines_reverse, "-", -(self.__genome_size - self.__offset_forward - 1))]

            for strand_content, strand, offset in prep:
                for position, c in strand_content.items():
                    nmeth, ntot = c["nmeth"], c["ntot"]
                    ratio = 0 if ntot == 0 else nmeth/ntot
                    csvfile.writerow([self.__chromoId, abs(offset + position), strand, c["context"], nmeth, ntot, ratio])
