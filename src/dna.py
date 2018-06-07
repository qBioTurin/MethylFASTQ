#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
from io import StringIO
import random
import enum

class strand(enum.Enum):
    Unspecified = 0
    ForwardStrand = 1
    ReverseStrand = 2
    BisulfiteForwardStrand = 3
    BisulfiteReverseStrand = 4
    BisulfiteForwardReverseComplement = 5
    BisulfiteReverseReverseComplement = 6

    #aggiungere metodi per fare il reverse, il complement etc e restituire il valore corretto



class dna(object):
    """Una sequenza di DNA è una sequenza di basi nucleotidiche"""

    def __init__(self, initial, begin=None, end=None):
        self.sequence = list(initial)
        self.__begin = begin
        self.__end = end

    @property
    def begin(self):
        return self.__begin

    @property
    def end(self):
        return self.__end

    def reverse_complement(self):
        return dna(map(self.__complement, self.sequence[::-1]), self.begin, self.end)

    def complement(self):
        return dna(map(self.__complement, self.sequence), self.begin, self.end)

    def reverse(self):
        return dna(self.sequence[::-1], self.begin, self.end)

    def bisulfite(self):
        return dna(map(self.__bisulfite, self.sequence), self.begin, self.end)

    def single_end_sequencing(self, read_length):
        return read(str(self)[:read_length])
#        return str(self)[:read_length]

    def paired_end_sequencing(self, read_length):
        r1 = dna(self[:read_length])
        r2 = dna(self[-read_length:]).reverse_complement()
        return paired_end_read(str(r1), str(r2))
#        return str(r1), str(r2)

    def __complement(self, base):
        if base == 'c':
            ret = 'g'
        elif base == 'C':
            ret = 'G'
        elif base == 'g':
            ret = 'c'
        elif base == 'G':
            ret = 'C'
        else:
            ret = 'a' if base.lower() == 't' else 't'

        return ret

    def __bisulfite(self, base):
        return base if base != 'C' else 't'

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return "".join(self.sequence)

    def __getitem__(self, key):
        return self.sequence[key]

    def __setitem__(self, key, value):
        self.sequence[key] = value

    def __delitem__(self, key):
        del self.sequence[key]


# class single_end_read(object):
#     """Una read single-end è un frammento single-stranded di dna,
#     a cui ad ogni base è associato un valore di qualità"""
#
#     def __init__(self, fragment, length):
#     #    dna.__init__(fragment.sequence)
#         self.__fragment = fragment.single_end_read(length)
#         self.__quality = [42] * len(self.__fragment)
#         self.__strand = fragment.strand
#
#     @property
#     def quality(self):
#         return self.__quality
#
#     @property
#     def strand(self):
#         return self.__strand
#
#     def set_sequencing_errors(self, error_rate):
#         """ """
#         for index in range(len(self.sequence)):
#             if random.uniform(0, 1) < error_rate:
#                 self.__sequence[index] = random.sample("actg", 1)[0]
#                 self.__quality[index] += 0
#         return self
#
#     def fastqize(self):
#         return "eheh"


class fragment(dna):
    """Un oggetto fragment è un frammento double-stranded di dna"""
    def __init__(self, sequence, begin, end):
        dna.__init__(self, sequence, begin, end)

    @property
    def forward_strand(self):
        return single_strand_fragment(self, strand.ForwardStrand)

    @property
    def reverse_strand(self):
        return single_strand_fragment(self.reverse_complement(), strand.ReverseStrand)

    def methylate(self, methyl_function): #parametri...
        """ fa cose e modifica la sequenza mettendo maiuscole o minuscole le c e le g... """

        sequence = [methyl_function(nt, pos) if nt in "cCGg" else nt for pos, nt in enumerate(self.sequence, self.begin)]
        return fragment(sequence, self.begin, self.end)


class single_strand_fragment(dna):
    def __init__(self, fragment, strand=strand.Unspecified):
        dna.__init__(self, fragment, fragment.begin, fragment.end)
        self.__strand = strand

    @property
    def strand(self):
        return self.__strand

    def bisulfite(self):
        new_strand = strand.BisulfiteForwardStrand if self.strand is strand.ForwardStrand else strand.BisulfiteReverseStrand
        return single_strand_fragment(fragment.bisulfite(self), strand=new_strand)

    def reverse_complement(self):
        if self.strand in (strand.ForwardStrand, strand.ReverseStrand):
            new_strand = strand.ReverseStrand if self.strand is strand.ForwardStrand else strand.ForwardStrand
        else:
            new_strand = strand.BisulfiteForwardReverseComplement if self.strand is strand.BisulfiteForwardStrand else strand.BisulfiteReverseReverseComplement
        dna_rc = super(single_strand_fragment, self).reverse_complement()
        return single_strand_fragment(dna_rc, strand=new_strand)

    def single_end_sequencing(self, read_length):
        return read(str(self)[:read_length], self.begin, self.strand)

    def paired_end_sequencing(self, read_length):
        r1 = read(str(self)[:read_length], self.begin, self.strand)
        r2 = dna(str(self)[-read_length:]).reverse_complement()
        r2 = read(str(r2), self.begin, self.strand)
        return paired_end_read(r1, r2)
#        r1 = dna(self[:read_length])41407932
#        r2 = dna(self[-read_length:]).reverse_complement()
#        return paired_end_read(str(r1), str(r2))


#    def single_end_read(self, length):
#        self.single_end_sequencing(length)
    #    return single_end_read(self, length)#.set_sequencing_errors(error_rate)
    #    return read2(self.).set_sequencing_errors(error_rate)
    # def paired_end_read(self, length):
    #     return paired_end_read(self, length)#.set_sequencing_errors(error_rate)

class read(object):
    #sequence è un cazzo di oggetto single_strand_fragment
    def __init__(self, sequence, begin, strand=strand.Unspecified):
        self.__sequence = dna(sequence, begin, 0)#list(sequence)
        self.quality = [42] * len(sequence) #default quality
        self.strand = strand #eheh
#        self.begin = 12

    @property
    def begin(self):
        return self.__sequence.begin

    def set_sequencing_errors(self, error_rate):
        """ """
        for index in range(len(self.sequence)):
            if random.uniform(0, 1) < error_rate:
                self.__sequence[index] = random.sample("actg", 1)[0]
                self.quality[index] -= 2
        return self

    def fastqize(self, read_id, fsize, offset):
        begin = self.begin + offset + 1 #

        #indici di mapping delle read
        #frammento forward:             [b, b+read_length]
        #frammento forward (undir):     uguale a frammento reverse
        #frammento reverse:             [b+fsize-rlength, b+fsize]
        #frammento reverse (undir):     uguale a frammento forward

        flags = "r"
        if self.strand in (strand.ForwardStrand, strand.BisulfiteForwardStrand, strand.BisulfiteForwardReverseComplement):
            flags = "f" #maremma doggy
        if self.strand in (strand.BisulfiteForwardReverseComplement, strand.BisulfiteReverseReverseComplement):
            flags += ":rc"
#        flags = "f" if self.strand in (strand.ForwardStrand, strand.BisulfiteForwardStrand, strand.BisulfiteForwardReverseComplement) else "r"

#        flags = "f" if self.strand is strand.ForwardStrand else "r"
#        if self.strand in (strand.ReverseStrand, strand.BisulfiteReverseStrand, strand.BisulfiteReverseReverseComplement):
        if self.strand in (strand.ReverseStrand, strand.BisulfiteReverseStrand, strand.BisulfiteForwardReverseComplement):
            begin += (fsize - len(self))
            
        read_id = "{}:{}:{}".format(read_id, begin, flags)
        default_quality = "~" * len(self)#self.__seqparams.read_length

        fastq_read = "@{}\n{}\n+\n{}\n".format(read_id, str(self.__sequence).upper(), default_quality)
        record = SeqIO.read(StringIO(fastq_read), "fastq-illumina")
        #set read quality
        record.letter_annotations["phred_quality"] = self.quality
        return record #Bio.SeqRecord.SeqRecord

    @property
    def sequence(self):
        return "".join(self.__sequence)

    def __len__(self):
        return len(self.__sequence)

class paired_end_read(object):
    def __init__(self, r1, r2):
        self.__r1 = r1#read(r1)
        self.__r2 = r2#read(r2)

    def set_sequencing_errors(self, error_rate):
        self.__r1.set_sequencing_errors(error_rate)
        self.__r2.set_sequencing_errors(error_rate)
        return self

    @property
    def r1(self):
        return self.__r1

    @property
    def r2(self):
        return self.__r2

    def fastqize(self, read_id, fsize, offset):
        #il nome delle read paired-end finisce con /1 o /2 -- modificare record fastq??
        r1 = self.__r1.fastqize(read_id, fsize, offset)
        r2 = self.__r2.fastqize(read_id, fsize, offset)

        return r1, r2
