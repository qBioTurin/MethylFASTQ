#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, argparse
import csv

def print_confusion_matrix(logfile, tool_name, cm):
    tp, fn, fp, tn = cm[0][0], cm[0][1], cm[1][0], cm[1][1]
    acc = (tp+tn) / (tp + fn + fp + tn)
    prec = tp / (tp + fp)
    rec = tp / (tp + fn)
    f1_score = 2 * prec * rec / (prec + rec)

    logfile.write("\nTool: {}\n".format(tool_name))
    logfile.write("TP: {}, FN: {}, FP: {}, TN: {}\n".format(tp, fn, fp, tn))
    logfile.write("Accuracy: {}\n".format(acc))
    logfile.write("Precision: {}\n".format(prec))
    logfile.write("Recall: {}\n".format(rec))
    logfile.write("F1 score: {}\n\n".format(f1_score))

def compare(gold_map, tool_map):
    cm = [[0,0],[0,0]]
    tp = set()

    for cytosine in gold_map.union(tool_map):
        gold_flag = cytosine in gold_map
        tool_flag = cytosine in tool_map
        cm[int(not gold_flag)][int(not tool_flag)] += 1
        if gold_flag and tool_flag:
            tp.add(cytosine)

    return tp, cm

def load_gold_file(filename):
    ext = filename.split(".")[-1].lower()
    cpg = set()

    with open(filename) as f:
        csvfile = csv.reader(f, delimiter="\t")

        if ext == "bed":
            for row in csvfile:
                ch, pos = row[:2]
                cytosine = ch, int(pos)
                cpg.add(cytosine)
        elif ext == "ch3":
            for row in csvfile:
                if row[3] == "CytosineContext.CG":
                    cytosine = (row[0], int(row[1]) + 1)
                    cpg.add(cytosine)
        else:
            sys.exit("Unexpected input file")

    return cpg

def load_bsmap_methylation(filename):
    cpg = set()

    with open(filename) as f:
        csvfile = csv.reader(f, delimiter="\t")
        next(csvfile, None) #skip header

        for row in csvfile:
            ch, pos, strand, context = row[:4]

            if context == "CG":
                cytosine = ch, int(pos)
                cpg.add(cytosine)

    return cpg

def load_bismark_methylation(filename):
    cpg = set()

    with open(filename) as f:
        for row in csv.reader(f, delimiter="\t"):
            ch, pos = row[2:4]
            cytosine = ch, int(pos)
            cpg.add(cytosine)

    return cpg


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="methyl-fastq")

    parser.add_argument("--gold", dest="gold_file", metavar="cg_gold_file", action="store", required=True,
                        type=str, help=".bed for real data or .ch3 for synthetic data")
    # parser.add_argument("--ch3", dest="ch3_file", metavar="ch_gold_file.ch3", action="store";
    #                     type=str, help="")
    parser.add_argument("--bsmap", dest="bsmap_file", metavar="methratio_output_file", required=True,
                        action="store", type=str, help="")
    parser.add_argument("--bismark", dest="bismark_file", metavar="bismark_methylation_extractor_output_file", required=True,
                        action="store", type=str, help="")
    parser.add_argument("--log", dest="logfile", metavar="log file", action="store", required=True, type=str, help="")

    args = parser.parse_args()

    logfile = open(args.logfile, "w")

    logfile.write("Reading CpG file {}...\n".format(args.gold_file))
    cpg_gold = load_gold_file(args.gold_file)


    # #legge file .ch3
    # with open(gold) as handle:
    #     csvfile = csv.reader(handle, delimiter="\t")
    #     for row in csvfile:
    #         if row[3] == "CytosineContext.CG":
    #             cytosine = (row[0], int(row[1]) + 1)
    #             cpg_map.add(cytosine)

    logfile.write("Analysing methratio file {}...\n".format(args.bsmap_file))
    bsmap_map = load_bsmap_methylation(args.bsmap_file)

    #legge output methratio
    # with open(methratio_output) as handle:
    #     csvfile = csv.reader(handle, delimiter="\t")
    #     next(csvfile, None) #skip header
    #
    #     for row in csvfile:
    #         ch, pos, strand, context = row[:4]
    #
    #         if context == "CG":
    #             cytosine = ch, int(pos)
    #             bsmap_map.add(cytosine)

    tp_bsmap, bsmap_confusion_matrix = compare(cpg_map, bsmap_map)
    print_confusion_matrix(logfile, "bsmap", bsmap_confusion_matrix)

    logfile.write("Analysing bismark methylation file {}...\n".format(args.bismark_file))
    bismark_map = load_bismark_methylation(args.bismark_file)


    #legge output bismark.methylation_extractor - file .Cpg_report.txt
    # with open(bismark_meth) as handle:
    #     csvfile = csv.reader(handle, delimiter="\t")
    # #    next(csvfile, None)
    #
    #     for row in csvfile:
    #         ch, pos = row[: 2]
    #         cytosine = ch, int(pos)
    #         bismark_map.add(cytosine)

    tp_bismark, bismark_confusion_matrix = compare(cpg_map, bismark_map)
    print_confusion_matrix(logfile, "bismark", bismark_confusion_matrix)

    in_common = len(tp_bismark.intersection(tp_bsmap))
    logfile.write("True Positives in common between BSMAP and Bismark: {}\n".format(in_common))
    logfile.close()
