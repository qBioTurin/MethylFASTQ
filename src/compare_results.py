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


if __name__ == "__main__":
    if len(sys.argv) < 5:
        sys.exit("Usage: {} cg_gold_file methratio_file bismark_file logfile".format(sys.argv[0]))

    gold = sys.argv[1]
    methratio_output = sys.argv[2]
    bismark_meth = sys.argv[3]
    logfile = open(sys.argv[4], "w")

    cpg_map = set()
    bsmap_map = set()
    bismark_map = set()

    logfile.write("Reading CpG file {}...\n".format(gold))

    #legge file .ch3
    with open(gold) as handle:
        csvfile = csv.reader(handle, delimiter="\t")
        for row in csvfile:
            if row[3] == "CytosineContext.CG":
                cytosine = (row[0], int(row[1]) + 1)
                cpg_map.add(cytosine)

    logfile.write("Analysing methratio file {}...\n".format(methratio_output))

    #legge output methratio
    with open(methratio_output) as handle:
        csvfile = csv.reader(handle, delimiter="\t")
        next(csvfile, None) #skip header

        for row in csvfile:
            ch, pos, strand, context = row[:4]

            if context == "CG":
                cytosine = ch, int(pos)
                bsmap_map.add(cytosine)

    tp_bsmap, bsmap_confusion_matrix = compare(cpg_map, bsmap_map)
    print_confusion_matrix(logfile, "bsmap", bsmap_confusion_matrix)

    logfile.write("Analysing bismark methylation file {}...\n".format(bismark_meth))


    #legge output bismark.methylation_extractor - file .Cpg_report.txt
    with open(bismark_meth) as handle:
        csvfile = csv.reader(handle, delimiter="\t")
    #    next(csvfile, None)

        for row in csvfile:
            ch, pos = row[: 2]
            cytosine = ch, int(pos)
            bismark_map.add(cytosine)

    tp_bismark, bismark_confusion_matrix = compare(cpg_map, bismark_map)
    print_confusion_matrix(logfile, "bismark", bismark_confusion_matrix)

    in_common = len(tp_bismark.intersection(tp_bsmap))
    logfile.write("True Positives in common between BSMAP and Bismark: {}\n".format(in_common))
    logfile.close()
