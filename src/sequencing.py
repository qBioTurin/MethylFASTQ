#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
from io import StringIO
import random, csv
import os, pickle
import multiprocessing as mp
from timeit import default_timer as timer
import enum
import sys

import dna


def format_time(time):
    strtime = "{:03.5f} m".format(time/60) if time >= 60 else "{:03.5f} s".format(time)
    return strtime

class CytosineContext(enum.Enum):
    CG = "CG"
    CHG = "CHG"
    CHH = "CHH"

    def __str__(self):
        return self.value

class Cytosine(object):
    def __init__(self, strand, context=CytosineContext.CHH):
        self.__context = context
        self.__strand = strand
        self.__nmeth = 0
        self.__ncov = 0

    @property
    def context(self):
        return self.__context

    @property
    def strand(self):
        return self.__strand

    @property
    def nmeth(self):
        return self.__nmeth

    @property
    def ncov(self):
        return self.__ncov

    def methylate(self):
        self.__nmeth += 1

    def covered(self):
        self.__ncov += 1


class Stats(object):
    def __init__(self):
        self.__num_reads = 0
        self.__covered_cs = 0
        self.__bp_sequenced = 0

    @property
    def nreads(self):
        return self.__num_reads

    @property
    def ncytosines(self):
        return self.__covered_cs

    @property
    def nbases(self):
        return self.__bp_sequenced

    def increment_reads(self, nr):
        self.__num_reads += nr

    def increment_cytosines(self, nc):
        self.__covered_cs += nc

    def increment_bps(self, nbp):
        self.__bp_sequenced += nbp

    def update(self, stats):
        self.increment_reads(stats.nreads)
        self.increment_cytosines(stats.ncytosines)
        self.increment_bps(stats.nbases)


class ChromosomeSequencer(object):
    """ ChromosomeSequencer legge un cromosoma da un file FASTA e
    scrive su un file FASTQ le read prodotte dal sequenziamento.
    È possibile scegliere il tipo di library (single-end o paired-end),
    se la library è direzionale o non direzionale, la lunghezza dei frammenti
    e delle read """

    def __init__(self, chromosome, target_regions=list()):
        """Parserizza il cromosoma individuando gli intervalli da sequenziare e scartando
        i nucleotidi indefiniti (basi N)"""

        chromosome_sequence = str(chromosome.seq).lower()

        self.__chromoId = chromosome.id
#        self.__genome_size = len(chromosome.seq)
        self.__fragments = [(chromosome_sequence[begin:end], begin, end) for (begin, end) in target_regions]
        self.__stats = Stats()

        if len(target_regions) == 0:
            #WGBS option
            print("Parsing {} sequence... ".format(chromosome.id), end="", flush=True)
        #    chromosome = str(chromosome.seq).lower()
            chr_size = len(chromosome.seq)

            #get salient fragments
            start, i = timer(), 0

            while i < chr_size:
                while i < chr_size and chromosome_sequence[i] == 'n':
                    i += 1
                begin = i
                while i < chr_size and chromosome_sequence[i] != 'n':
                    i += 1

                t = (chromosome_sequence[begin:i], begin, i)
                self.__fragments.append(t)
            else:
                last = self.__fragments.pop()
                if last[1] < last[2]:
                    self.__fragments.append(last)

            tot_time = format_time(timer() - start)

        #sommo dimensioni intervalli
        self.__stats.increment_bps(sum([e-b for _, b, e in self.__fragments]))
        print("{} fragments found: {} bp. Elapsed time {}".format(len(self.__fragments), self.__stats.nbases, tot_time), flush=True)


    def sequencing(self, params):
        queue = mp.Queue() #comunicazione tra processi figli e padre

        #ordino frammenti in modo decrescente rispetto alla loro lunghezza
        self.__fragments = [{"seq": x[0], "from": x[1], "to": x[2], "params": params} for x in self.__fragments]
        self.__fragments.sort(key=lambda el: el["to"] - el["from"], reverse=True)

        print([x["to"] - x["from"] for x in self.__fragments])

        single_end = params.seq_mode == "single_end"

        #inizializzo i processi con i relativi input: parametri vari (frammento e offset) + coda + indice per il join
        processes = [mp.Process(target=self.create_reads, name="methylPIPPO-{}".format(index), args=(param, params.seq_mode, (index, queue))) \
                        for index, param in enumerate(self.__fragments)]


        num_input = len(self.__fragments)   #numero di input da processare
        num_process = params.num_processes  #numero massimo di processi da utilizzare

        curr = 0        #indice del prossimo processo da startare
        curr_exec = 0   #numero processi attualmente in esecuzione

        #starto primi processi
        while curr_exec < min(num_input, num_process):
            processes[curr].start()
            curr += 1
            curr_exec += 1


        output_filename = self.__get_output_filename(params)

        if single_end:
            fastq_file = open("{}.fastq".format(output_filename), "a")
        else:
            fastq_file1 = open("{}_R1.fastq".format(output_filename), "a")
            fastq_file2 = open("{}_R2.fastq".format(output_filename), "a")

        meth_file = open("{}.ch3".format(output_filename), "a")
        csv_meth = csv.writer(meth_file, delimiter="\t")

        #finchè tutti i processi non sono stati eseguiti e non hanno terminato...
        while not (curr == num_input and curr_exec == 0):
            val = queue.get()
            tval = type(val)

            if tval is int: #segnale di terminazione di un processo
                processes[val].join()

                if curr < num_input:
                    processes[curr].start()
                    curr += 1
                else:
                    curr_exec -= 1
            elif tval is tuple:
                tag, data = val
                dsize = len(data)

                if tag == "fastq_se":
                    for record in data:
                        SeqIO.write(record, fastq_file, "fastq")
                    else:
                        self.__stats.increment_reads(dsize)

                elif tag == "fastq_pe":
                    for read1, read2 in data:
                        SeqIO.write(read1, fastq_file1, "fastq")
                        SeqIO.write(read2, fastq_file2, "fastq")
                    else:
                        self.__stats.increment_reads(dsize)

                elif tag == "ch3":
                    for record in data:
                        csv_meth.writerow(record)
                    else:
                        self.__stats.increment_cytosines(dsize)
        else: #fine while -> chiudo file
            meth_file.close()

            if single_end:
                fastq_file.close()
            else:
                fastq_file1.close()
                fastq_file2.close()

        for p in processes:
            p.join()

        return self.__stats

    def __get_output_filename(self, params):
        fasta = "".join(params.fasta_file.split("/")[-1].split(".")[:-1])
        se_pe = "se" if params.seq_mode == "single_end" else "pe"
        dir_nondir = "dir" if params.lib_mode == "directional" else "undir"
        return "{}/{}_{}_f{}r{}_{}".format(params.output_path, fasta, se_pe, params.fragment_size, params.read_length, dir_nondir)


    def create_reads(self, data, seq_mode, index_queue):
        offset_begin, offset_end = data["from"], data["to"]
        sequence = data["seq"]
        params = data["params"]
        index, queue = index_queue

        print("processing fragment from {} to {}".format(data["from"], data["to"]))
        fs = FragmentSequencer(self.__chromoId, sequence, offset_begin, offset_end, \
                                params, queue=queue)

        fs.single_end_sequencing() if seq_mode == "single_end" else fs.paired_end_sequencing()

        queue.put(index)



class FragmentSequencer(object):
    """..."""

    def __init__(self, chr_id, sequence, begin_seq, end_seq, seqparams, queue):
        #informazioni sulla sequenza
        self.__chromoId = chr_id
        self.__sequence = sequence
#        self.__genome_size = gsize
        #posizione del frammento nel genoma
        self.__offset_forward = begin_seq
#        self.__offset_reverse = gsize - end_seq
        #dimensione buffer + parametri vari
        self.__buffer_size = seqparams.buffer_size
        self.__seqparams = seqparams
        self.__p_meth =  {
            CytosineContext.CG: seqparams.p_cg,
            CytosineContext.CHG: seqparams.p_chg,
            CytosineContext.CHH: seqparams.p_chh
        }
        #dati sulle citosine
        self.__cytosines = dict()
        self.__queue = queue

        self.__read_quality = dna.read_quality(seqparams.read_length, seqparams.min_quality, seqparams.max_quality)
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
            fastq_fw = bsfs.single_end_sequencing(read_length, self.__read_quality)\
                        .set_sequencing_errors(error_rate)\
                        .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

            bsrs = fragment.reverse_strand.bisulfite()
            fastq_rv = bsrs.single_end_sequencing(read_length, self.__read_quality)\
                        .set_sequencing_errors(error_rate)\
                        .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

            reads.extend([fastq_fw, fastq_rv])

            if not directional:
                fastq_fwrc = bsfs.reverse_complement()\
                            .single_end_sequencing(read_length, self.__read_quality)\
                            .set_sequencing_errors(error_rate)\
                            .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

                fastq_rvrc = bsrs.reverse_complement()\
                            .single_end_sequencing(read_length, self.__read_quality)\
                            .set_sequencing_errors(error_rate)\
                            .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

                reads.extend([fastq_fwrc, fastq_rvrc])


            if len(reads) > self.__buffer_size:
                self.__queue.put(("fastq_se", reads.copy()))
    #            self.__persist_record(reads)
                reads.clear()
        else:
            if len(reads) > 0:
                self.__queue.put(("fastq_se", reads.copy()))
    #            self.__persist_record(reads)
                reads.clear()

        cinfo = self.__format_methylation()
        self.__queue.put(("ch3", cinfo))

#        self.__persist_methylation()


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
            fastq_fw = bsfs.paired_end_sequencing(read_length, self.__read_quality)\
                        .set_sequencing_errors(error_rate)\
                        .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

            bsrs = fragment.reverse_strand.bisulfite()
            fastq_rv = bsrs.paired_end_sequencing(read_length, self.__read_quality)\
                        .set_sequencing_errors(error_rate)\
                        .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

            reads.extend([fastq_fw, fastq_rv])



            if not directional:
                fastq_fwrc = bsfs.reverse_complement()\
                            .paired_end_sequencing(read_length, self.__read_quality)\
                            .set_sequencing_errors(error_rate)\
                            .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

                fastq_rvrc = bsrs.reverse_complement()\
                            .paired_end_sequencing(read_length, self.__read_quality)\
                            .set_sequencing_errors(error_rate)\
                            .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

                reads.extend([fastq_fwrc, fastq_rvrc])

            if len(reads) > self.__buffer_size:
                self.__queue.put(("fastq_pe", reads.copy()))
        #        self.__persist_record(reads)
                reads.clear()
        else:
            if len(reads) > 0:
                self.__queue.put(("fastq_pe", reads.copy()))
        #        self.__persist_record(reads)
                reads.clear()

        cinfo = self.__format_methylation()
        self.__queue.put(("ch3", cinfo))
    #    self.__persist_methylation()

    ###################### Methylation ######################

    def __initialize_cytosines(self):
        """Parserizza il genoma e indicizza le citosine sui due strand"""

        limit = len(self.__sequence)

        for pos, base in enumerate(self.__sequence):
            #strand +
            if base == 'c':
                context = CytosineContext.CHH
                if pos+1 < limit and self.__sequence[pos+1] == 'g':
                    context = CytosineContext.CG
                elif pos+2 < limit and self.__sequence[pos+2] == 'g':
                    context = CytosineContext.CHG

                self.__cytosines[pos] = Cytosine("+", context)
            #strand -
            elif base == 'g':
                context = CytosineContext.CHH
                if pos-1 >= 0 and self.__sequence[pos-1] == "c":
                    context = CytosineContext.CG
                elif pos-2 >= 0 and self.__sequence[pos-2] == "c":
                    context = CytosineContext.CHG

                self.__cytosines[pos] = Cytosine("-", context)


    def methylate_cytosine(self, base, position):
        state = base.upper()
        if position in self.__cytosines:
            cytosine = self.__cytosines[position]
            cytosine.covered()

            if random.uniform(0, 1) <= self.__p_meth[cytosine.context]:
                state = state.lower()
                cytosine.methylate()

        return state


    ###################### Introduzione mutazioni ######################

    def __set_snp(self):
        """Setta SNP random sulla reference"""

        self.__sequence = "".join([random.sample("actg", 1)[0] if random.uniform(0, 1) < self.__seqparams.snp_rate else base\
                                   for base in self.__sequence])


    def __format_methylation(self):
        """ Returns a list containing all the information about cytosines """
        beta_score = lambda c: 0 if c.ncov == 0 else c.nmeth / c.ncov

        return [(self.__chromoId, abs(self.__offset_forward+position), cytosine.strand, cytosine.context,\
                cytosine.nmeth, cytosine.ncov, beta_score(cytosine)) \
                    for position, cytosine in self.__cytosines.items()]
